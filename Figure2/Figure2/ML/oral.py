from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import RepeatedStratifiedKFold, train_test_split
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from imblearn.over_sampling import SMOTE
import seaborn as sns


def cross_validate_model(model, X_train, y_train):
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=42)

    fold_accuracies = []
    fold_precisions = []
    fold_recalls = []
    fold_f1_scores = []
    confusion_matrices = []

    for fold, (train_index, val_index) in enumerate(rskf.split(X_train, y_train)):
        X_train_fold, X_val_fold = X_train.iloc[train_index], X_train.iloc[val_index]
        y_train_fold, y_val_fold = y_train.iloc[train_index], y_train.iloc[val_index]

        model.fit(X_train_fold, y_train_fold)

        y_pred = model.predict(X_val_fold)

        accuracy = accuracy_score(y_val_fold, y_pred)
        precision = precision_score(y_val_fold, y_pred, average='weighted')
        recall = recall_score(y_val_fold, y_pred, average='weighted')
        f1 = f1_score(y_val_fold, y_pred, average='weighted')

        fold_accuracies.append(accuracy)
        fold_precisions.append(precision)
        fold_recalls.append(recall)
        fold_f1_scores.append(f1)

        cm = confusion_matrix(y_val_fold, y_pred)
        confusion_matrices.append(cm)

    mean_accuracy = np.mean(fold_accuracies)
    mean_precision = np.mean(fold_precisions)
    mean_recall = np.mean(fold_recalls)
    mean_f1_score = np.mean(fold_f1_scores)

    summed_cm = np.sum(confusion_matrices, axis=0)

    return {
        "accuracy": mean_accuracy,
        "precision": mean_precision,
        "recall": mean_recall,
        "f1_score": mean_f1_score,
        "confusion_matrix": summed_cm
    }


def evaluate_on_holdout(model, X_train, y_train, X_holdout, y_holdout):
    model.fit(X_train, y_train)

    y_pred = model.predict(X_holdout)

    accuracy = accuracy_score(y_holdout, y_pred)
    precision = precision_score(y_holdout, y_pred, average='weighted')
    recall = recall_score(y_holdout, y_pred, average='weighted')
    f1 = f1_score(y_holdout, y_pred, average='weighted')

    cm = confusion_matrix(y_holdout, y_pred)

    return {
        "accuracy": accuracy,
        "precision": precision,
        "recall": recall,
        "f1_score": f1,
        "confusion_matrix": cm
    }


data_folder = Path("data")
df = pd.read_csv(data_folder / "oral_machine.csv")
X = df.iloc[:, 2:]
y = df["Group"]
ros = SMOTE(random_state=20)
X_resampled, y_resampled = ros.fit_resample(X, y)
X_train, X_holdout, y_train, y_holdout = train_test_split(X_resampled, y_resampled, test_size=0.4, stratify=y_resampled, random_state=42)
classifiers = {
    "Random Forest": RandomForestClassifier(random_state=20, n_estimators=200, class_weight="balanced", max_depth=5, min_samples_split=5),
    "SVM": SVC(random_state=20, C=9,probability=True, class_weight="balanced"),
    "Gradient Boosting": GradientBoostingClassifier(random_state=50, max_depth=6, subsample=0.4),
    "Logistic Regression": LogisticRegression(max_iter=500, random_state=55)
}
cross_val_results = {}
holdout_results = {}
holdout_predictions = {}
for name, model in classifiers.items():
    print(f"Cross-validating {name}...")
    cross_val_results[name] = cross_validate_model(model, X_train, y_train)

    print(f"Evaluating {name} on hold-out set...")
    holdout_results[name] = evaluate_on_holdout(model, X_train, y_train, X_holdout, y_holdout)

    y_holdout_prob = model.predict_proba(X_holdout)
    holdout_predictions[name] = y_holdout_prob
for name in classifiers.keys():
    print(f"\n{name} Results from Cross-Validation:")
    print(f"Accuracy: {cross_val_results[name]['accuracy']:.4f}")
    print(f"Precision: {cross_val_results[name]['precision']:.4f}")
    print(f"Recall: {cross_val_results[name]['recall']:.4f}")
    print(f"F1 Score: {cross_val_results[name]['f1_score']:.4f}")

    plt.figure(figsize=(8, 6))
    sns.heatmap(cross_val_results[name]["confusion_matrix"], annot=True, fmt='d', cmap='Blues',
                xticklabels=np.unique(y), yticklabels=np.unique(y))
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title(f'{name} Cross-Validation Confusion Matrix')
    plt.show()

    print(f"\n{name} Results on Hold-Out Set:")
    print(f"Accuracy: {holdout_results[name]['accuracy']:.4f}")
    print(f"Precision: {holdout_results[name]['precision']:.4f}")
    print(f"Recall: {holdout_results[name]['recall']:.4f}")
    print(f"F1 Score: {holdout_results[name]['f1_score']:.4f}")

    plt.figure(figsize=(8, 6))
    sns.heatmap(holdout_results[name]["confusion_matrix"], annot=True, fmt='d', cmap='Blues', xticklabels=np.unique(y),
                yticklabels=np.unique(y))
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title(f'{name} Hold-Out Set Confusion Matrix')
    plt.show()


holdout_results_df = pd.DataFrame(columns=[
    "Model", "Accuracy", "Precision", "Recall", "F1 Score"
])


for name, results in holdout_results.items():
    new_row = pd.DataFrame({
        "Model": [name],
        "Accuracy": [results["accuracy"]],
        "Precision": [results["precision"]],
        "Recall": [results["recall"]],
        "F1 Score": [results["f1_score"]]
    })
    holdout_results_df = pd.concat([holdout_results_df, new_row], ignore_index=True)

holdout_results_csv_path = "holdout_model_evaluation_oral_results.csv"
holdout_results_df.to_csv(holdout_results_csv_path, index=False)

print(f"Hold-Out  {holdout_results_csv_path}")
