

mkdir -p rabbitmq_data


chmod 777 -R rabbitmq_data


mkdir -p core/templates/
mkdir -p core/static/
mkdir -p core/uploads/

rm -rf core/*/migrations/0*.py
rm -rf core/*/migrations/1*.py
rm -rf core/*/migrations/2*.py
rm -rf core/*/migrations/3*.py
rm -rf core/*/migrations/4*.py
rm -rf core/*/migrations/5*.py
rm -rf core/*/migrations/6*.py
rm -rf core/*/migrations/7*.py
rm -rf core/*/migrations/8*.py
rm -rf core/*/migrations/9*.py

