FROM nvidia/cuda:12.6.3-cudnn-devel-ubuntu24.04


# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

ENV TZ=America/Chicago

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1


# Basic update, upgrade, and install some packages
RUN apt-get update -y
RUN apt-get upgrade -y

RUN apt-get install -y --no-install-recommends build-essential

RUN apt-get install -y git wget curl gpg pkg-config default-libmysqlclient-dev \
libxml2-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
libbz2-dev libffi-dev libsqlite3-dev liblzma-dev libreadline-dev


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

RUN python3 -m pip install Django psycopg2-binary mysqlclient sqlalchemy python-dotenv graphviz whitenoise gunicorn celery pika django-celery-results djangorestframework django-cors-headers

RUN python3 -m pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

RUN python3 -m pip install requests pandas scanpy anndata scipy performer_pytorch scikit-learn


# Changed working directory
WORKDIR /home/scplant_backend/core/
