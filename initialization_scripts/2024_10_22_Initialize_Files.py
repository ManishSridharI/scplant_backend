import os
import requests

# Base URL of the API
BASE_URL = "http://digbio-g2pdeep.rnet.missouri.edu:8449"
# BASE_URL = "http://127.0.0.1:8449"

# Login credentials
LOGIN_URL = f"{BASE_URL}/accounts/api/login/"
USERNAME = "chanye"
PASSWORD = "Qwerty#0987"

# API endpoint to make calls
ORGANISM_API_ENDPOINT = f"{BASE_URL}/organisms/api/organism_generate/"
PREDICTOR_API_ENDPOINT = f"{BASE_URL}/predictors/api/predictor_generate/"
H5AD_DATASET_API_ENDPOINT = f"{BASE_URL}/h5addatasets/api/h5ad_dataset_generate/"
SCRIPT_API_ENDPOINT = f"{BASE_URL}/scripts/api/script_generate/"

BASE_DIR = "/data/html/Prod/scPlant/scplant_backend/"
UPLOAD_DIR = BASE_DIR + "core/uploads/"

# Organism list
ORGANISM_ARRAY = ["Arabidopsis Thaliana", "Zea Mays", "Oryza Sativa", "Glycine Max"]

# Predictor list
PREDICTOR_DICT = {
    "arabidopsis_hvg20k": "model_codebase/models/arabidopsis_hvg20k.ckpt",
    "corn_hvg10k": "model_codebase/models/corn_hvg10k.ckpt",
    "rice_hvg20k": "model_codebase/models/rice_hvg20k.ckpt",
    "soybean_hvg20k": "model_codebase/models/soybean_hvg20k.ckpt"
}

PREDICTOR_ORGANISM_DICT = {
    "arabidopsis_hvg20k": 1,
    "corn_hvg10k": 2,
    "rice_hvg20k": 3,
    "soybean_hvg20k": 4
}

# Dataset list
H5AD_DATASET_DICT = {
    "corn_SRP335_hvg10k": "model_codebase/public_data/corn_SRP335_hvg10k.h5ad",
    "rice_SRP286_hvg20k": "model_codebase/public_data/rice_SRP286_hvg20k.h5ad",
    "soybean_flowerbud_hvg20k": "model_codebase/public_data/soybean_flowerbud_hvg20k.h5ad",
    "SRP171_hvg20k": "model_codebase/public_data/SRP171_hvg20k.h5ad",
    "SRP235_hvg20k": "model_codebase/public_data/SRP235_hvg20k.h5ad",
    "SRP330_hvg20k": "model_codebase/public_data/SRP330_hvg20k.h5ad",
    "tester": "model_codebase/user_datasets/tester.h5ad"
}

H5AD_DATASET_ORGANISM_DICT = {
    "corn_SRP335_hvg10k": 2,
    "rice_SRP286_hvg20k": 3,
    "soybean_flowerbud_hvg20k": 4,
    "SRP171_hvg20k": 1,
    "SRP235_hvg20k": 1,
    "SRP330_hvg20k": 1,
    "tester": 1
}

# Python script list
PYTHON_SCRIPT_DICT = {
    "inference": "model_codebase/inference.py",
    "annotate_and_plot": "model_codebase/annotate_and_plot.py",
    "control_vs_treatment": "model_codebase/control_vs_treatment.py",
    "compare_celltype_distributions": "model_codebase/compare_celltype_distributions.py"
}

# R script list
R_SCRIPT_DICT = {
    "convert_rds_to_h5ad": "model_codebase/convert_rds_to_h5ad.R"
}

# Start a session to persist cookies (e.g., CSRF token)
session = requests.Session()

# Step 1: Login and get CSRF token
login_data = {
    "username": USERNAME,
    "password": PASSWORD
}

# Perform login
login_response = session.post(LOGIN_URL, data=login_data)

# Check for successful login
if login_response.status_code == 200:
    print("Login successful!")
else:
    print("Login failed!", login_response.text)
    exit()

# Extract CSRF token (depends on how the server sends it)
csrf_token = session.cookies.get("csrftoken")  # Assuming token is in cookies
if not csrf_token:
    print("CSRF token not found!")
    exit()

# Step 2: Make API call with CSRF token
headers = {
    "X-CSRFToken": csrf_token,  # Include CSRF token in header
    "Content-Type": "application/json"  # JSON content type
}

for element in ORGANISM_ARRAY:
    payload = {
        "organism_name": element
    }
    response = session.post(
        ORGANISM_API_ENDPOINT,
        json=payload,
        headers=headers
    )
    if response.status_code == 201:
        print("API call successful!")
        print("Response:", response.json())
    else:
        print("API call failed!")
        print("Status Code:", response.status_code)
        print("Response:", response.text)

predictor_output_array = []
h5ad_dataset_output_array = []
python_script_output_array = []
r_script_output_array = []

for key in PREDICTOR_DICT.keys():
    payload = {
        "predictor_name": key,
        "predictor_filename": key,
        "predictor_organism": PREDICTOR_ORGANISM_DICT[key],
        "predictor_public_flag": True
    }
    response = session.post(
        PREDICTOR_API_ENDPOINT,
        json=payload,
        headers=headers
    )
    if response.status_code == 201:
        predictor_output_array.append(
            "sudo cp -rf " + BASE_DIR + PREDICTOR_DICT[key] + " " + UPLOAD_DIR + str(response.json()['Predictor'])
        )
        print("API call successful!")
        print("Response:", response.json())
    else:
        print("API call failed!")
        print("Status Code:", response.status_code)
        print("Response:", response.text)

for key in H5AD_DATASET_DICT.keys():
    payload = {
        "h5ad_dataset_name": key,
        "h5ad_dataset_file_extension": "h5ad",
        "h5ad_dataset_filename": key,
        "h5ad_dataset_organism": H5AD_DATASET_ORGANISM_DICT[key],
        "h5ad_dataset_public_flag": True
    }
    response = session.post(
        H5AD_DATASET_API_ENDPOINT,
        json=payload,
        headers=headers
    )
    if response.status_code == 201:
        h5ad_dataset_output_array.append(
            "sudo cp -rf " + BASE_DIR + H5AD_DATASET_DICT[key] + " " + UPLOAD_DIR + str(response.json()['H5adDataset'])
        )
        print("API call successful!")
        print("Response:", response.json())
    else:
        print("API call failed!")
        print("Status Code:", response.status_code)
        print("Response:", response.text)

for key in PYTHON_SCRIPT_DICT.keys():
    payload = {
        "script_name": key,
        "script_filename": key,
        "script_file_extension": "py",
        "script_public_flag": True
    }
    response = session.post(
        SCRIPT_API_ENDPOINT,
        json=payload,
        headers=headers
    )
    if response.status_code == 201:
        python_script_output_array.append(
            "sudo cp -rf " + BASE_DIR + PYTHON_SCRIPT_DICT[key] + " " + UPLOAD_DIR + str(response.json()['Script'])
        )
        if key == "inference" or key == "annotate_and_plot":
            python_script_output_array.append(
                "sudo cp -rf " + BASE_DIR + "model_codebase/performer_pytorch" + " " + UPLOAD_DIR + os.path.dirname(
                    str(response.json()['Script'])) + "/"
            )
            python_script_output_array.append(
                "sudo cp -rf " + BASE_DIR + "model_codebase/utils.py" + " " + UPLOAD_DIR + os.path.dirname(
                    str(response.json()['Script'])) + "/"
            )
        print("API call successful!")
        print("Response:", response.json())
    else:
        print("API call failed!")
        print("Status Code:", response.status_code)
        print("Response:", response.text)

for key in R_SCRIPT_DICT.keys():
    payload = {
        "script_name": key,
        "script_filename": key,
        "script_file_extension": "R",
        "script_public_flag": True
    }
    response = session.post(
        SCRIPT_API_ENDPOINT,
        json=payload,
        headers=headers
    )
    if response.status_code == 201:
        r_script_output_array.append(
            "sudo cp -rf " + BASE_DIR + R_SCRIPT_DICT[key] + " " + UPLOAD_DIR + str(response.json()['Script'])
        )
        print("API call successful!")
        print("Response:", response.json())
    else:
        print("API call failed!")
        print("Status Code:", response.status_code)
        print("Response:", response.text)

with open(str(BASE_DIR + "initialization_scripts/2024_10_23_Copy_Files.sh"), "w") as writer:
    writer.write("\n\n")
    writer.write("sudo chmod 775 -R " + str(UPLOAD_DIR))
    writer.write("\n\n")
    writer.write("\n".join(predictor_output_array))
    writer.write("\n\n")
    writer.write("\n".join(h5ad_dataset_output_array))
    writer.write("\n\n")
    writer.write("\n".join(python_script_output_array))
    writer.write("\n\n")
    writer.write("\n".join(r_script_output_array))
    writer.write("\n\n")
    writer.write("sudo chmod 775 -R " + str(UPLOAD_DIR))
    writer.write("\n\n")
