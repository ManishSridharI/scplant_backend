
base_dir="/data/html/Prod/scPlant/scplant_backend"

sudo chmod 775 -R ${base_dir}/core/uploads/*

mkdir -p ${base_dir}/core/uploads/predictors/chanye/
mkdir -p ${base_dir}/core/uploads/datasets/chanye/
mkdir -p ${base_dir}/core/uploads/scripts/chanye/

cp -rf ${base_dir}/model_codebase/models/* ${base_dir}/core/uploads/predictors/chanye/
cp -rf ${base_dir}/model_codebase/public_data/* ${base_dir}/core/uploads/datasets/chanye/
cp -rf ${base_dir}/model_codebase/user_datasets/* ${base_dir}/core/uploads/datasets/chanye/
cp -rf ${base_dir}/model_codebase/*.py ${base_dir}/core/uploads/scripts/chanye/
cp -rf ${base_dir}/model_codebase/*.R ${base_dir}/core/uploads/scripts/chanye/

ls -lha core/uploads/*/chanye/
