from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from sklearn import svm
from sklearn.model_selection import KFold, StratifiedKFold
from utils import *
from plot import *
import joblib


# 加载数据
def load_data(data_folder, csv_name):
    # 加载特征数据，跳过第一行表头和第一列
    filename = os.path.join(data_folder, csv_name)
    cell_name = []
    features = []

    with open(filename, 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        # Skip header
        header = next(csv_reader)

        for row in csv_reader:
            cell_name.append(row[0])
            features.append(list(map(float, row[1:])))

    return cell_name, np.array(features)


def load_labels(label_file_path):
    # 加载标签数据，使用第二列作为标签
    labels = np.loadtxt(label_file_path, delimiter=',', skiprows=1, usecols=2)

    return labels


# 训练和验证模型，一次十折交叉实验
def train_and_validate(cell_name, features, labels, model_name, exp_num, random_seed):

    if model_name == "Lasso":
        model = linear_model.Lasso()
    elif model_name == "LinearSVR":
        model = svm.SVR(kernel='linear')
    elif model_name == "LinearSVR":
        model = svm.SVR()
    elif model_name == "BayesianRidge":
        model = linear_model.BayesianRidge()
    elif model_name == "RandomForest":
        model = RandomForestRegressor(n_estimators=30, random_state=0)
    else:
        print("Model name error")
        return
        # model = linear_model.LinearRegression()

    MAEs = []
    MSEs = []
    RMSEs = []
    r2s = []
    all_train_preds = []
    all_test_preds = []

    print("features: ", features.shape)
    print("labels: ", labels.shape)

    label_file_path = os.path.join(data_folder, "labels.csv")

    # skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=random_seed)
    skf = KFold(n_splits=10, shuffle=True, random_state=random_seed)
    for n_fold, (train_idx, test_idx) in enumerate(skf.split(features, labels)):
        print(f"EXP{exp_num}: Fold{n_fold}, Training ... ")
        X_train, y_train = features[train_idx], labels[train_idx]
        X_test, y_test = features[test_idx], labels[test_idx]

        model.fit(X_train, y_train)
        y_train_pred = model.predict(X_train)

        output_folder = os.path.join(f'../checkpoint/HC', f'{model_name}')
        os.makedirs(output_folder, exist_ok=True)

        train_output_filename = os.path.join(output_folder, f'train_score_exp{exp_num}_fold{n_fold}.csv')
        write_to_csv(train_output_filename, np.array(cell_name)[train_idx], y_train_pred, y_train)
        plot_regression(train_output_filename, output_folder, label_file_path)

        joblib.dump(model, os.path.join(output_folder, f'{model_name}_exp{exp_num}_fold{n_fold}.model'))

        print("Validating ... ")
        y_test_pred = model.predict(X_test)

        all_train_preds.append(y_train_pred)
        all_test_preds.append(y_test_pred)

        output_filename = os.path.join(output_folder, f'test_score_exp{exp_num}_fold{n_fold}.csv')
        write_to_csv(output_filename, np.array(cell_name)[test_idx], y_test_pred, y_test)
        plot_regression(output_filename, output_folder, label_file_path)

        mae = mean_absolute_error(y_test, y_test_pred)
        mse = mean_squared_error(y_test, y_test_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_test_pred))
        r2 = r2_score(y_test, y_test_pred)

        MAEs.append(mae)
        MSEs.append(mse)
        RMSEs.append(rmse)
        r2s.append(r2)

        print(f"### EXP{exp_num} ### Fold{n_fold} ### MAE={mae} ### MSE={mse} ### RMSE={rmse} ### R2={r2}")

        # 独立测试集
        print("Testing ... ")
        test_output_folder = f'../checkpoint/Disease/{model_name}'
        os.makedirs(test_output_folder, exist_ok=True)
        test_cell_name, test_features = load_data(data_folder, 'test_data.csv')
        test_label_file_path = os.path.join(data_folder, "test_labels.csv")
        test_labels = load_labels(test_label_file_path)
        test_pred = model.predict(test_features)

        mae = mean_absolute_error(test_labels, test_pred)
        mse = mean_squared_error(test_labels, test_pred)
        rmse = np.sqrt(mean_squared_error(test_labels, test_pred))
        print(f"### EXP{exp_num} ### Fold{n_fold} ### MAE={mae} ### MSE={mse} ### RMSE={rmse}")

        output_filename = os.path.join(test_output_folder, f'score_exp{exp_num}_fold{n_fold}.csv')
        write_to_csv(output_filename, test_cell_name, test_pred, test_labels)

        plot_regression(output_filename, test_output_folder, test_label_file_path)

    return MAEs, MSEs, RMSEs, r2s, all_train_preds, all_test_preds


def main(data_folder, random_seeds, model_name_list):
    for i, seed in enumerate(random_seeds):
        print(f"===== 第{i}次实验，使用随机种子: {seed} 进行实验 =====")

        label_file_path = os.path.join(data_folder, "labels.csv")
        labels = load_labels(label_file_path)

        # 加载数据
        cell_name, features = load_data(data_folder, 'data.csv')

        # 评价指标
        MAEs_avg = []
        MSEs_avg = []
        RMSEs_avg = []
        r2s_avg = []

        # 模型
        for model_name in model_name_list:
            MAEs, MSEs, RMSEs, r2s, all_train_preds, all_test_preds = train_and_validate(cell_name, features, labels, model_name, i, seed)
            MAEs_avg.extend(MAEs)
            MSEs_avg.extend(MSEs)
            RMSEs_avg.extend(RMSEs)
            r2s_avg.extend(r2s)

            print(f"=====  {model_name} 第{i+1}次实验结果 =====")
            print("MAE: {:.4}".format(np.mean(MAEs)), "R2: {:.4}".format(np.mean(r2s)))
            print("MSE: {:.4}".format(np.mean(MSEs)), "RMSE: {:.4}".format(np.mean(RMSEs)))

            with open(f"../checkpoint/{model_name}_{i}_{seed}.txt", "w") as f:
                print("========= 实验结果 =========", file=f)
                print("MAE {:.2f} ± {:.2f}".format(np.mean(MAEs_avg), np.std(MAEs_avg)), file=f)
                print("MSE {:.2f} ± {:.2f}".format(np.mean(MSEs_avg), np.std(MSEs_avg)), file=f)
                print("RMSE {:.2f} ± {:.2f}".format(np.mean(RMSEs_avg), np.std(RMSEs_avg)), file=f)
                print("R2 {:.2f} ± {:.2f}".format(np.mean(r2s_avg), np.std(r2s_avg)), file=f)


if __name__ == "__main__":
    data_folder = f'../data'
    result_folder = '../checkpoint'
    os.makedirs(result_folder, exist_ok=True)

    num_experiments = 10
    file_path = "../checkpoint/seeds.txt"
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            random_seeds = [int(seed) for seed in f.read().split()]
    else:
        random_seeds = np.random.randint(0, 500, num_experiments)
        with open(file_path, "w") as f:
            print(*random_seeds, file=f)

    # model_name_list = ["Lasso", "RandomForest", "LinearSVR", "SVR", "BayesianRidge"]
    model_name_list = ["RandomForest"]
    print(random_seeds)
    print(model_name_list)
    main(data_folder, random_seeds, model_name_list)
