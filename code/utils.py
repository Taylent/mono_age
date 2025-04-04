from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import csv
import math


# 计算评价指标并保存到文件
def calculate_and_save_metrics(labels, preds, experiment_num, feature_name):
    mae = mean_absolute_error(labels, preds)
    r2 = r2_score(labels, preds)
    rmse = math.sqrt(mean_squared_error(labels, preds))

    metrics_file_path = f'../checkpoint/{feature_name}_exp{experiment_num}.txt'
    with open(metrics_file_path, 'w') as f:
        f.write(f"MAE: {mae}\n")
        f.write(f"R2: {r2}\n")
        f.write(f"RMSE: {rmse}\n")


def write_to_csv(filename, subject_id, preds, labels):
    data = []
    for index, pred in enumerate(preds):
        sample_name = subject_id[index]
        label = str(labels[index])
        pred = str(pred)
        data.append([sample_name, label, pred])

    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['cell_name', 'true_value', 'predicted_value'])
        writer.writerows(data)