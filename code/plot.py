import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from itertools import cycle
import os
from sklearn.linear_model import LinearRegression


def plot_regression(csv_file, png_save_file, label_file_path):
    # Load labels.csv and result.csv
    labels_df = pd.read_csv(label_file_path)
    results_df = pd.read_csv(csv_file)

    # Merge dataframes on 'sample_name'
    merged_df = pd.merge(results_df, labels_df, on='cell_name')

    # Group by 'group' and plot regression curves
    groups = merged_df['disease_type'].unique()
    colors = cycle(plt.cm.tab10.colors)  # Cycle through colors for different groups

    base_file_name = os.path.splitext(os.path.basename(csv_file))[0]

    plt.figure(figsize=(12, 6))

    for group in groups:
        group_data = merged_df[merged_df['disease_type'] == group]

        X = group_data['true_value'].values
        Y = group_data['predicted_value'].values

        # Fit a linear regression model
        model = LinearRegression()
        X_train = group_data['true_value'].values.reshape(-1, 1)
        model.fit(X_train, Y)

        # 将线性模型的拟合结果画成虚线
        color = next(colors)
        plt.plot(X_train, model.predict(X_train), color=color, linestyle='-', linewidth=2, label=f'Regression Model ({group})')
        # 将模型预测的结果画成散点
        plt.scatter(X, Y, color=color, label=f'{group}')

    # 生成x数据
    x = np.linspace(0, 100, 800)
    y = x
    plt.plot(x, y, 'r--', label='y=x')

    plt.xlabel('True Age')
    plt.ylabel('Pre Age')
    plt.title(f'Regression Curves for singleCell')
    plt.legend()
    plt.savefig(os.path.join(png_save_file, base_file_name+"_age.png"))
    plt.close()

    for group in groups:
        group_data = merged_df[merged_df['disease_type'] == group]

        X = group_data['true_value'].values
        Y = group_data['predicted_value'].values

        # Plot true_age vs pre_age - true_age for the current group
        plt.scatter(X, Y-X, color=next(colors), label=f'{group}')

    plt.xlabel('True Age')
    plt.ylabel('Scores')
    plt.title(f'Regression Curves for singleCell')
    plt.legend()
    plt.savefig(os.path.join(png_save_file, base_file_name+"_score.png"))
    plt.close()


if __name__ == '__main__':
    data_folder = f'../data'
    model_name = "RandomForest"
    for exp_num in range(10):
        for n_fold in range(10):
            test_output_folder = f'../checkpoint/Disease/{model_name}'
            output_filename = os.path.join(test_output_folder, f'score_exp{exp_num}_fold{n_fold}_test.csv')
            test_label_file_path = os.path.join(data_folder, "test_labels.csv")
            plot_regression(output_filename, test_output_folder, test_label_file_path)