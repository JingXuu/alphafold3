CUDA_VISIBLE_DEVICES=2 bash complex/code/run_af3_fold_add_binder_split.sh \
    ./complex/fold_try \
    ./data/split_msas/af3_datapipeline_output/af3-dp-20251221-021152-272270354 \
    ~/public_databases/models \
    ./data/af3_fold_results

CUDA_VISIBLE_DEVICES=2 nohup bash /home/xujing/alphafold3/complex/code/run_af3_fold_add_binder_split.sh \
    /home/xujing/BindCraft/results/for_mpnn_3.5 \
    ./data/split_msas/af3_datapipeline_output/af3-dp-20251221-021152-272270354 \
    ~/public_databases/models \
    ./data/mpnn35_results \
    > ./data/logs/MPNN35/af3_6KL5_7DY7_8EWV.log 2>&1 &

# > ./data/logs/MPNN35/af3_6KL5_7DY7_8EWV_$(date +%Y%m%d-%H%M%S).log 2>&1 &