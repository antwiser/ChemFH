#!/bin/bash
hyperopt_dir=checkpoint/DMPNN/hyperopt
model_dir=checkpoint/DMPNN/$1
train_path=FH_modeldata/model_data/$1/train.csv
val_path=FH_modeldata/model_data/$1/valid.csv
test_path=FH_modeldata/model_data/$1/test.csv
dataset_type=classification
metric=auc
gpu=1

# Hyperparameter optimization
chemprop_hyperopt \
    --dataset_type $dataset_type \
    --data_path $train_path \
    --separate_val_path $val_path \
    --separate_test_path $val_path \
    --gpu $gpu \
    --batch_size 128 \
    --num_iters 30 \
    --epochs 300 \
    --metric $metric \
    --class_balance \
    --num_workers 40 \
    --aggregation norm \
    --search_parameter_keywords depth ffn_num_layers  hidden_size ffn_hidden_size dropout \
    --config_save_path $hyperopt_dir/config.json \
    --hyperopt_checkpoint_dir $hyperopt_dir \
    --log_dir $hyperopt_dir


# Training with optimized hyperparameters
chemprop_train \
    --dataset_type $dataset_type \
    --data_path $train_path \
    --separate_val_path $val_path \
    --separate_test_path $test_path \
    --epochs 300 \
    --gpu $gpu \
    --class_balance \
    --config_path $hyperopt_dir/config.json \
    --save_dir $model_dir \
    --ensemble_size 1 \
    --save_preds \
    --extra_metrics accuracy \
    --show_individual_scores 
