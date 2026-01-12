import pandas as pd
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import joblib
import os

# Configuration
CSV_PATH = "../cache_results.csv"
MODEL_SAVE_PATH = "fitness_model.pth"
SCALER_X_PATH = "scaler_x.joblib"
SCALER_Y_PATH = "scaler_y.joblib"

def load_data(csv_path):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"File not found: {csv_path}")
    
    df = pd.read_csv(csv_path)
    
    # Input columns (10 params)
    # Based on the CSV header: gamma_d_GB_GluSynapse,gamma_p_GB_GluSynapse,a00,a01,a10,a11,a20,a21,a30,a31
    feature_cols = [
        'gamma_d_GB_GluSynapse', 'gamma_p_GB_GluSynapse',
        'a00', 'a01', 'a10', 'a11', 'a20', 'a21', 'a30', 'a31'
    ]
    
    # Target columns (6 scores)
    all_target_cols = [
        'fitness_score_0', 'fitness_score_1', 'fitness_score_2', 
        'fitness_score_3', 'fitness_score_4', 'fitness_score_5'
    ]
    
    # Check which target columns have data
    valid_target_cols = []
    for col in all_target_cols:
        if col in df.columns and df[col].notna().sum() > 0:
            valid_target_cols.append(col)
        else:
            print(f"Warning: Column {col} is entirely NaN or missing. Skipping.")
            
    if not valid_target_cols:
        raise ValueError("No valid target columns found!")
        
    # Drop rows with NaNs in feature cols or valid target cols
    df = df.dropna(subset=feature_cols + valid_target_cols)
    
    X = df[feature_cols].values
    y = df[valid_target_cols].values
    
    return X, y, feature_cols, valid_target_cols

class FitnessPredictor(nn.Module):
    def __init__(self, input_dim, output_dim):
        super(FitnessPredictor, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 64),
            nn.ReLU(),
            nn.Linear(64, output_dim)
        )
        
    def forward(self, x):
        return self.net(x)

def train_model():
    print(f"Loading data from {CSV_PATH}...")
    try:
        X, y, feature_cols, target_cols = load_data(CSV_PATH)
    except FileNotFoundError:
        # Fallback for running from different CWD
        alt_path = os.path.join(os.path.dirname(__file__), "../cache_results.csv")
        print(f"Trying alternative path: {alt_path}")
        X, y, feature_cols, target_cols = load_data(alt_path)

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    # Normalize
    scaler_x = StandardScaler()
    X_train_scaled = scaler_x.fit_transform(X_train)
    X_test_scaled = scaler_x.transform(X_test)
    
    scaler_y = StandardScaler()
    y_train_scaled = scaler_y.fit_transform(y_train)
    y_test_scaled = scaler_y.transform(y_test)
    
    # Convert to tensors
    X_train_tensor = torch.FloatTensor(X_train_scaled)
    y_train_tensor = torch.FloatTensor(y_train_scaled)
    X_test_tensor = torch.FloatTensor(X_test_scaled)
    y_test_tensor = torch.FloatTensor(y_test_scaled)
    
    # Initialize model
    model = FitnessPredictor(input_dim=len(feature_cols), output_dim=len(target_cols))
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    
    # Training loop
    epochs = 2000
    print("Starting training...")
    best_val_loss = float('inf')
    
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad()
        outputs = model(X_train_tensor)
        loss = criterion(outputs, y_train_tensor)
        loss.backward()
        optimizer.step()
        
        if (epoch + 1) % 100 == 0:
            model.eval()
            with torch.no_grad():
                val_outputs = model(X_test_tensor)
                val_loss = criterion(val_outputs, y_test_tensor)
            
            print(f'Epoch [{epoch+1}/{epochs}], Loss: {loss.item():.4f}, Val Loss: {val_loss.item():.4f}')
            
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                torch.save(model.state_dict(), MODEL_SAVE_PATH)
            
    # Save scalers and metadata
    joblib.dump(scaler_x, SCALER_X_PATH)
    joblib.dump(scaler_y, SCALER_Y_PATH)
    joblib.dump(target_cols, "target_cols.joblib")
    print(f"Best model saved to {MODEL_SAVE_PATH}")
    print(f"Scalers saved to {SCALER_X_PATH} and {SCALER_Y_PATH}")
    print(f"Target columns saved to target_cols.joblib")

if __name__ == "__main__":
    train_model()
