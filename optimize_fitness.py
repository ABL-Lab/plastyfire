import torch
import numpy as np
import joblib
import pandas as pd
from scipy.optimize import minimize
import os
import sys

# Import from local script
try:
    from train_fitness_model import FitnessPredictor, CSV_PATH, MODEL_SAVE_PATH, SCALER_X_PATH, SCALER_Y_PATH, load_data
except ImportError:
    # If running from a different directory, we might need to adjust path
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from train_fitness_model import FitnessPredictor, CSV_PATH, MODEL_SAVE_PATH, SCALER_X_PATH, SCALER_Y_PATH, load_data

def load_model_and_scalers():
    if not os.path.exists(MODEL_SAVE_PATH):
        raise FileNotFoundError(f"Model file not found: {MODEL_SAVE_PATH}. Run train_fitness_model.py first.")
        
    # Load scalers
    scaler_x = joblib.load(SCALER_X_PATH)
    scaler_y = joblib.load(SCALER_Y_PATH)
    
    # Load target columns
    if os.path.exists("target_cols.joblib"):
        target_cols = joblib.load("target_cols.joblib")
    else:
        # Fallback
        target_cols = [f'fitness_score_{i}' for i in range(6)]
        
    # Load model
    input_dim = 10
    output_dim = len(target_cols)
    
    model = FitnessPredictor(input_dim, output_dim)
    model.load_state_dict(torch.load(MODEL_SAVE_PATH))
    model.eval()
    
    return model, scaler_x, scaler_y, target_cols

def objective_function(params, model, scaler_x, scaler_y):
    # params is a 1D numpy array of shape (10,)
    
    # Scale input
    params_reshaped = params.reshape(1, -1)
    params_scaled = scaler_x.transform(params_reshaped)
    params_tensor = torch.FloatTensor(params_scaled)
    
    # Predict
    with torch.no_grad():
        predicted_scaled = model(params_tensor)
    
    # Inverse scale output
    predicted_scores = scaler_y.inverse_transform(predicted_scaled.numpy())
    
    # Objective: Sum of fitness scores
    # We want to minimize this sum
    total_score = np.sum(predicted_scores)
    return total_score

def optimize():
    print("Loading model...")
    model, scaler_x, scaler_y, target_cols = load_model_and_scalers()
    
    # Initial guess (mean of data or random)
    print(f"Loading data for initial guess from {CSV_PATH}...")
    try:
        df = pd.read_csv(CSV_PATH)
    except FileNotFoundError:
        alt_path = os.path.join(os.path.dirname(__file__), "../cache_results.csv")
        df = pd.read_csv(alt_path)

    feature_cols = [
        'gamma_d_GB_GluSynapse', 'gamma_p_GB_GluSynapse',
        'a00', 'a01', 'a10', 'a11', 'a20', 'a21', 'a30', 'a31'
    ]
    
    # Use the parameters that gave the minimum sum so far as starting point
    # Only use available target cols
    valid_cols = [c for c in target_cols if c in df.columns]
    
    # Drop NaNs
    df = df.dropna(subset=feature_cols + valid_cols)
    
    df['sum_score'] = df[valid_cols].sum(axis=1)
    best_idx = df['sum_score'].idxmin()
    initial_guess = df.loc[best_idx, feature_cols].values
    best_existing_sum = df.loc[best_idx, 'sum_score']
    
    print(f"Initial guess (best from data): {initial_guess}")
    print(f"Initial sum of scores (from data): {best_existing_sum}")
    print(f"Optimizing for sum of: {target_cols}")
    
    # Bounds
    bounds = []
    for col in feature_cols:
        min_val = df[col].min()
        max_val = df[col].max()
        # Allow some extrapolation (e.g., 20% beyond observed range)
        range_val = max_val - min_val
        margin = range_val * 0.2
        bounds.append((min_val - margin, max_val + margin))
    
    print("Starting optimization (L-BFGS-B)...")
    result = minimize(
        objective_function, 
        initial_guess, 
        args=(model, scaler_x, scaler_y),
        method='L-BFGS-B',
        bounds=bounds
    )
    
    print("\nOptimization Result:")
    print(f"Success: {result.success}")
    print(f"Message: {result.message}")
    print(f"Optimal Parameters: {result.x}")
    print(f"Predicted Minimum Sum of Scores: {result.fun}")
    
    if result.fun < best_existing_sum:
        print(f"IMPROVEMENT: Found parameters predicted to be better than best data point by {best_existing_sum - result.fun:.4f}")
    else:
        print("No improvement found over best data point (model might be conservative or local minimum).")
    
    # Verify with model prediction breakdown
    params_reshaped = result.x.reshape(1, -1)
    params_scaled = scaler_x.transform(params_reshaped)
    with torch.no_grad():
        pred_scaled = model(torch.FloatTensor(params_scaled))
    pred_scores = scaler_y.inverse_transform(pred_scaled.numpy())[0]
    
    print("\nPredicted Scores Breakdown:")
    for i, score in enumerate(pred_scores):
        col_name = target_cols[i] if i < len(target_cols) else f"score_{i}"
        print(f"{col_name}: {score:.4f}")
        
    # Save result to file
    result_df = pd.DataFrame([result.x], columns=feature_cols)
    result_df['predicted_sum_score'] = result.fun
    for i, score in enumerate(pred_scores):
        col_name = target_cols[i] if i < len(target_cols) else f"score_{i}"
        result_df[f'predicted_{col_name}'] = score
        
    result_df.to_csv("optimized_params.csv", index=False)
    print("\nSaved optimized parameters to optimized_params.csv")

if __name__ == "__main__":
    optimize()
