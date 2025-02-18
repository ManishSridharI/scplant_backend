# FROM ubuntu
# FROM nvidia/cuda:12.6.3-cudnn-devel-ubuntu24.04
FROM nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04


# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

ENV TZ=America/Chicago

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1


# Basic update, upgrade, and install some packages
RUN apt-get update -y
RUN apt-get upgrade -y

RUN apt-get install -y --no-install-recommends build-essential

RUN apt-get install -y lsb-release software-properties-common

RUN apt-get install -y git wget curl gpg pkg-config default-libmysqlclient-dev \
libxml2-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
libbz2-dev libffi-dev libsqlite3-dev liblzma-dev libreadline-dev \
libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
automake


# Setup Pyenv and Python then install Python packages
RUN curl https://pyenv.run | bash

RUN echo '' >> ~/.bashrc
RUN echo '' >> ~/.bashrc
RUN echo 'export PYENV_ROOT="/root/.pyenv"' >> ~/.bashrc
RUN echo 'command -v pyenv >/dev/null || export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
RUN echo 'eval "$(pyenv init -)"' >> ~/.bashrc
RUN echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.bashrc

ENV PYENV_ROOT="/root/.pyenv"
ENV PATH="$PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH"

RUN pyenv update
RUN pyenv install 3.12
RUN pyenv global 3.12

RUN python3 -m pip install --upgrade pip

RUN python3 -m pip install Django psycopg2-binary mysqlclient sqlalchemy \
python-dotenv graphviz whitenoise gunicorn celery pika django-celery-results \
djangorestframework django-cors-headers

RUN python3 -m pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

RUN python3 -m pip install requests pandas numpy scanpy anndata scipy \
performer_pytorch scikit-learn xlsxwriter openpyxl plotly

# RUN python3 -c "import torch; print(torch.__version__); print(torch.cuda.is_available())"


# Setup R and R packages
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

RUN gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

RUN add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

RUN apt-get -y install r-base

RUN R -e "install.packages('packrat', dependencies=TRUE,  repos='https://cloud.r-project.org/')"


# Changed working directory
WORKDIR /home/scplant_backend/core/
