#!/bin/bash
hyperopt_dir=checkpoint/DMPNN_Des/hyperopt
model_dir=checkpoint/DMPNN_Des/$1

train_path=FH_modeldata/model_data/$1/train.csv
val_path=FH_modeldata/model_data/$1/valid.csv
test_path=FH_modeldata/model_data/$1/test.csv

features_path=FH_modeldata/Des/$1/train_2d.npy
separate_val_features_path=FH_modeldata/Des/$1/valid_2d.npy
separate_test_features_path=FH_modeldata/Des/$1/test_2d.npy

out_test=output/test_Des.csv 
out_train=output/train_Des.csv 
out_valid=output/valid_Des.csv 

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
    --features_path $features_path \
    --separate_val_features_path $separate_val_features_path \
    --separate_test_features_path $separate_test_features_path \
    --no_features_scaling \
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
    --dataset_type classification \
    --data_path $train_path \
    --separate_val_path $val_path \
    --separate_test_path $test_path \
    --epochs 300\
    --gpu $gpu \
    --class_balance \
    --config_path $hyperopt_dir/config.json \
    --save_dir $model_dir \
    --ensemble_size 1 \
    --save_preds \
    --extra_metrics accuracy \
    --show_individual_scores \
    --no_features_scaling \
    --aggregation norm \
    --features_path $features_path \
    --separate_val_features_path $separate_val_features_path \
    --separate_test_features_path $separate_test_features_path


