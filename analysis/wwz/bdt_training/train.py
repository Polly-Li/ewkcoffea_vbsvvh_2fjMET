import matplotlib.pyplot as plt
import numpy as np
import mpltern
import pickle
import gzip
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
import xgboost as xgb

# Avoid import but not used error for mpltern:
#     - We need to import it in order to register the "ternary" projection
#     - But mpltern is not explicitly used anywhere, making flake8 upset
#     - Thus just add an assert line, so that flake8 will not catch it
assert mpltern

dd = pickle.load(gzip.open("bdt.pkl.gz"))

X_train_of = dd["X_train_of"]
y_train_of = dd["y_train_of"]
w_train_of = dd["w_train_of"]
X_train_sf = dd["X_train_sf"]
y_train_sf = dd["y_train_sf"]
w_train_sf = dd["w_train_sf"]
X_test_of = dd["X_test_of"]
y_test_of = dd["y_test_of"]
w_test_of = dd["w_test_of"]
X_test_sf = dd["X_test_sf"]
y_test_sf = dd["y_test_sf"]
w_test_sf = dd["w_test_sf"]
y_train_tmva_of = dd["y_train_tmva_of"]
y_test_tmva_of = dd["y_test_tmva_of"]
y_train_tmva_sf = dd["y_train_tmva_sf"]
y_test_tmva_sf = dd["y_test_tmva_sf"]

xgb_clf = xgb.XGBClassifier(objective='multi:softmax',
                            num_class=3,
                            missing=1,
                            booster="gbtree",
                            grow_policy="depthwise",
                            learning_rate=2.0,
                            n_estimators=100,
                            eval_metric=['merror','mlogloss'],
                            device="cuda",
                            seed=42,
                            n_jobs=32)

params = xgb_clf.get_params()
print("XGBoost Classifier Parameters:\n")
for param, value in params.items():
    print(f"{param}: {value}")

xgb_clf.fit(X_train_of,
            y_train_of,
            sample_weight=np.maximum(w_train_of, 0),
            verbose=1, # set to 1 to see xgb training round intermediate results
            eval_set=[(X_train_of, y_train_of), (X_test_of, y_test_of)])

xgb_clf.save_model("bdt.json")

# preparing evaluation metric plots
results = xgb_clf.evals_result()
epochs = len(results['validation_0']['mlogloss'])
x_axis = range(0, epochs)

# xgboost 'mlogloss' plot
fig, ax = plt.subplots(figsize=(9,5))
ax.plot(x_axis, results['validation_0']['mlogloss'], label='Train')
ax.plot(x_axis, results['validation_1']['mlogloss'], label='Test')
ax.legend()
plt.ylabel('mlogloss')
plt.title('GridSearchCV XGBoost mlogloss')
plt.savefig("mlogloss.pdf")

# xgboost 'merror' plot
fig, ax = plt.subplots(figsize=(9,5))
ax.plot(x_axis, results['validation_0']['merror'], label='Train')
ax.plot(x_axis, results['validation_1']['merror'], label='Test')
ax.legend()
plt.ylabel('merror')
plt.title('GridSearchCV XGBoost merror')
plt.savefig("merror.pdf")

y_prob_of = xgb_clf.predict_proba(X_test_of)
n_classes = len(np.unique(y_test_of))
fpr = {}
tpr = {}
roc_auc = {}
y_test_of_bin = label_binarize(y_test_of, classes=np.unique(y_test_of))
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test_of_bin[:, i], y_prob_of[:, i], sample_weight=np.abs(w_test_of))
    roc_auc[i] = auc(fpr[i], tpr[i])

# Plot ROC curve for each class
plt.figure(figsize=(8, 6))
for i in range(n_classes):
    plt.plot(fpr[i], tpr[i], lw=2, label='ROC curve (class %d) (AUC = %0.2f)' % (i, roc_auc[i]))

plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.savefig("roc.pdf")

# Ternary plot
fig = plt.figure(figsize=(25, 25))
ax0 = fig.add_subplot(2, 2, 1, projection="ternary")
ax1 = fig.add_subplot(2, 2, 2, projection="ternary")
ax2 = fig.add_subplot(2, 2, 3, projection="ternary")
ax3 = fig.add_subplot(2, 2, 4, projection="ternary")

ax0.tribin(y_prob_of[:, 0][y_test_of == 0], y_prob_of[:, 1][y_test_of == 0], y_prob_of[:, 2][y_test_of == 0], gridsize=20 , edgecolors="none", mincnt=-1)
ax1.tribin(y_prob_of[:, 0][y_test_of == 1], y_prob_of[:, 1][y_test_of == 1], y_prob_of[:, 2][y_test_of == 1], gridsize=20 , edgecolors="none", mincnt=-1)
ax2.tribin(y_prob_of[:, 0][y_test_of == 2], y_prob_of[:, 1][y_test_of == 2], y_prob_of[:, 2][y_test_of == 2], gridsize=20 , edgecolors="none", mincnt=-1)

ax3.scatter(y_prob_of[:, 0][y_test_of == 0], y_prob_of[:, 1][y_test_of == 0], y_prob_of[:, 2][y_test_of == 0], s=1 , marker='o' , color=(1 , 0 , 0 , 0.55))
ax3.scatter(y_prob_of[:, 0][y_test_of == 1], y_prob_of[:, 1][y_test_of == 1], y_prob_of[:, 2][y_test_of == 1], s=1 , marker='o' , color=(0 , 0 , 1 , 0.55))
ax3.scatter(y_prob_of[:, 0][y_test_of == 2], y_prob_of[:, 1][y_test_of == 2], y_prob_of[:, 2][y_test_of == 2], s=1 , marker='o' , color=(0 , 1 , 0 , 0.55))
plt.savefig("scatter.png")


#### TMVA

y_prob_of = y_test_tmva_of
n_classes = len(np.unique(y_test_of))
fpr = {}
tpr = {}
roc_auc = {}
y_test_of_bin = label_binarize(y_test_of, classes=np.unique(y_test_of))
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test_of_bin[:, i], y_prob_of[:, i], sample_weight=np.abs(w_test_of))
    roc_auc[i] = auc(fpr[i], tpr[i])

# Plot ROC curve for each class
plt.figure(figsize=(8, 6))
for i in range(n_classes):
    plt.plot(fpr[i], tpr[i], lw=2, label='ROC curve (class %d) (AUC = %0.2f)' % (i, roc_auc[i]))

plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.savefig("roctmva.pdf")

