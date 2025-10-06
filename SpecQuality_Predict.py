import pandas as pd
import numpy as np
import joblib
import os

# === CLASS DEFINITIONS ===
class EnsembleClassifier:
    def __init__(self, models, scalers, feature_columns):
        self.models = models
        self.scalers = scalers
        self.feature_columns = feature_columns
        
    def predict_proba(self, X):
        # Ensure we use the same features and order
        X = X[self.feature_columns]
        
        all_predictions = []
        for model, scaler in zip(self.models, self.scalers):
            X_scaled = scaler.transform(X)
            pred = model.predict_proba(X_scaled)[:, 1]  # Probability for class 1
            all_predictions.append(pred)
        
        # Average predictions across all models
        ensemble_proba = np.mean(all_predictions, axis=0)
        return np.column_stack([1 - ensemble_proba, ensemble_proba])
    
    def predict(self, X, threshold=0.5):
        proba = self.predict_proba(X)[:, 1]
        return (proba >= threshold).astype(int)

class EnsembleRegressor:
    def __init__(self, models, scalers, feature_columns):
        self.models = models
        self.scalers = scalers
        self.feature_columns = feature_columns
        
    def predict(self, X):
        X = X[self.feature_columns]
        
        all_predictions = []
        for model, scaler in zip(self.models, self.scalers):
            X_scaled = scaler.transform(X)
            pred = model.predict(X_scaled)
            all_predictions.append(pred)
        
        return np.mean(all_predictions, axis=0)

def calculate_geometric_mean(features, feature_set):
    """
    Calculate geometric mean for specified feature set without normalization
    """
    # Select only the features in the current feature set
    available_features = [f for f in feature_set if f in features.columns]
    top_features = features[available_features]
    
    # Calculate geometric mean directly without normalization
    geom_mean = np.exp(np.mean(np.log(top_features + 1e-10), axis=1))
    
    return geom_mean

def load_model(model_path):
    """
    Load a single model with proper error handling
    """
    print(f"Loading model from {model_path}...")
    
    try:
        model = joblib.load(model_path)
        print(f"✅ Successfully loaded model: {os.path.basename(model_path)}")
        
        # Determine model type and extract components
        if hasattr(model, 'predict') and hasattr(model, 'feature_columns'):
            # It's an ensemble model
            model_type = 'ensemble'
            feature_columns = model.feature_columns
            scaler = None
            is_classifier = hasattr(model, 'predict_proba')
            
        elif isinstance(model, dict):
            # It's a components dictionary
            model_type = 'components'
            feature_columns = model.get('feature_columns', [])
            scaler = model.get('scaler')
            is_classifier = 'clf_model' in model
            
            # Extract the actual model
            if is_classifier:
                model = model.get('clf_model')
            else:
                model = model.get('reg_model')
        else:
            print(f"❌ Unknown model type in {model_path}")
            return None, None, None, None
        
        print(f"   Model type: {'Classifier' if is_classifier else 'Regressor'}")
        print(f"   Features: {len(feature_columns)}")
        
        return model, feature_columns, scaler, is_classifier
        
    except Exception as e:
        print(f"❌ Error loading model {model_path}: {e}")
        return None, None, None, None

def predict_with_model(model, feature_columns, scaler, data, is_classifier):
    """
    Make predictions using a loaded model
    """
    # Check if all required features are available
    missing_features = [f for f in feature_columns if f not in data.columns]
    if missing_features:
        print(f"❌ Missing features: {missing_features}")
        return None, None
    
    # Extract features
    X = data[feature_columns]
    
    # Scale if scaler is available
    if scaler is not None:
        X_processed = scaler.transform(X)
    else:
        X_processed = X
    
    # Make predictions
    try:
        if is_classifier:
            probabilities = model.predict_proba(X_processed)[:, 1]
            predictions = (probabilities >= 0.5).astype(int)
            return predictions, probabilities
        else:
            predictions = model.predict(X_processed)
            return predictions, None
    except Exception as e:
        print(f"❌ Prediction error: {e}")
        return None, None

def calculate_sqs_scores(input_file, model_paths, output_dir=None):
    """
    Calculate SQS scores using multiple models on one input file
    """
    # Load input data
    print(f"Loading data from {input_file}...")
    df = pd.read_csv(input_file, sep='\t')
    print(f"✅ Loaded {len(df)} spectra")
    
    # Prepare output filename
    input_basename = os.path.basename(input_file)
    rootname = os.path.splitext(input_basename)[0]
    if output_dir is None:
        output_dir = os.path.dirname(input_file)
    output_file = os.path.join(output_dir, f"{rootname}_SQS.tsv")
    
    # Initialize results
    results = df.copy()
    
    # Calculate Geometric Mean scores (without class predictions)
    print("\nCalculating Geometric Mean scores...")
    
    # Define feature sets
    SQS5_FEATURES = ['Intensity', 'ComplementsFraction', 'PeakDensity', 'PeakCount', 'GoodDiffFraction']
    SQS10_FEATURES = ['PeakCount', 'SignalPeaksCount', 'SNR', 'Intensity', 'PeakDensity', 
                     'GoodDiffFraction', 'ComplementsFraction', 'IsotopePeaks', 
                     'NeutralLossPeaks', 'AverageRelativeIntensity']
    
    # Calculate geometric means
    results['SQS5_gm'] = calculate_geometric_mean(df, SQS5_FEATURES)
    results['SQS10_gm'] = calculate_geometric_mean(df, SQS10_FEATURES)
    
    # Process each model
    print("\nProcessing models...")
    model_results = {}
    
    for i, model_path in enumerate(model_paths):
        if not os.path.exists(model_path):
            print(f"❌ Model file not found: {model_path}")
            continue
            
        model_name = os.path.basename(model_path).replace('.pkl', '')
        print(f"\n🔧 Processing {model_name}...")
        
        # Load model
        model, feature_columns, scaler, is_classifier = load_model(model_path)
        if model is None:
            continue
        
        # Make predictions
        predictions, probabilities = predict_with_model(model, feature_columns, scaler, df, is_classifier)
        
        if predictions is not None:
            # Add results with model-specific column names
            if is_classifier:
                results[f'{model_name}_prob'] = probabilities
                results[f'{model_name}_class'] = predictions
                model_results[model_name] = {
                    'type': 'classifier',
                    'high_quality': predictions.sum(),
                    'percent_high_quality': predictions.sum() / len(predictions) * 100
                }
            else:
                results[f'{model_name}_score'] = predictions
                model_results[model_name] = {
                    'type': 'regressor',
                    'min_score': predictions.min(),
                    'max_score': predictions.max(),
                    'mean_score': predictions.mean()
                }
            print(f"✅ Added predictions from {model_name}")
        else:
            print(f"❌ Failed to get predictions from {model_name}")
    
    # Save results
    print(f"\nSaving results to {output_file}...")
    results.to_csv(output_file, sep='\t', index=False)
    
    # Print summary statistics
    print("\n" + "="*60)
    print("SQS CALCULATION SUMMARY")
    print("="*60)
    print(f"Total spectra processed: {len(results):,}")
    
    print(f"\nGeometric Mean Scores:")
    print(f"  SQS5 Geometric Mean:  {results['SQS5_gm'].min():.3f} to {results['SQS5_gm'].max():.3f}")
    print(f"  SQS10 Geometric Mean: {results['SQS10_gm'].min():.3f} to {results['SQS10_gm'].max():.3f}")
    
    if model_results:
        print(f"\nModel Predictions:")
        for model_name, stats in model_results.items():
            if stats['type'] == 'classifier':
                print(f"  {model_name}: {stats['high_quality']:,} high quality spectra ({stats['percent_high_quality']:.1f}%)")
            else:
                print(f"  {model_name}: {stats['min_score']:.3f} to {stats['max_score']:.3f} (mean: {stats['mean_score']:.3f})")
    
    # Show output columns
    #new_columns = [col for col in results.columns if col not in df.columns]
    #print(f"\nNew columns added: {new_columns}")
    
    return results

def batch_calculate_sqs(input_files, model_paths, output_dir):
    """
    Calculate SQS scores for multiple files
    """
    all_results = {}
    
    for input_file in input_files:
        if not os.path.exists(input_file):
            print(f"Warning: File {input_file} not found, skipping...")
            continue
            
        print(f"\n{'='*50}")
        print(f"Processing: {input_file}")
        print(f"{'='*50}")
        
        try:
            results = calculate_sqs_scores(input_file, model_paths, output_dir)
            if results is not None:
                all_results[input_file] = results
                print(f"✅ Successfully processed {input_file}")
            else:
                print(f"❌ Failed to process {input_file}")
                
        except Exception as e:
            print(f"❌ Error processing {input_file}: {str(e)}")
            import traceback
            traceback.print_exc()
    
    return all_results

# Example usage
if __name__ == "__main__":
    # Configuration - You can use 1-4 models here
    MODEL_PATHS = [  
        "SQS5_reg.pkl",       # Regression model
        "SQS5_clf.pkl",       # Classification model
        "SQS10_reg.pkl",      # Regression model
        "SQS10_clf.pkl",      # Classification model
    ]
    
    # Filter out models that don't exist (use at least one)
    available_models = [model for model in MODEL_PATHS if os.path.exists(model)]
    
    if not available_models:
        print("❌ No model files found! Please check the file paths.")
        exit(1)
    
    print(f"🎯 Using {len(available_models)} available models:")
    for model in available_models:
        print(f"   - {os.path.basename(model)}")
    
    # Input file
    #INPUT_FILE = "ad_pl01_filtered_new_contlabel.tsv"
    #INPUT_FILE = "NEW_QTOF_all_SQS_V7_MSGF_contlabel.tsv"
    INPUT_FILE = "UWA419_A_filtered_new_SQS_pyV7_contlabel.tsv"
    
    OUTPUT_DIR = "."
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Calculate SQS scores with multiple models
    results = calculate_sqs_scores(INPUT_FILE, available_models, OUTPUT_DIR)
    
    if results is not None:
        print(f"\n🎉 Success! Processed {len(results)} spectra with {len(available_models)} models.")
        #print("\nFirst 3 rows of results:")
        #new_cols = [col for col in results.columns if col not in ['SQS5_gm', 'SQS10_gm']]
        #print(results[new_cols].head(3))
    else:
        print("\n💥 Processing failed!")