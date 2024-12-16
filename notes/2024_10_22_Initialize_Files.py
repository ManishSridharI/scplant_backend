import requests


# Base URL of the API
BASE_URL = "http://digbio-g2pdeep.rnet.missouri.edu:8449"

# Login credentials
LOGIN_URL = f"{BASE_URL}/accounts/api/login/"
USERNAME = "chanye"
PASSWORD = "Qwerty#0987"

# API endpoint to make calls
PREDICTOR_API_ENDPOINT = f"{BASE_URL}/predictors/api/predictor_generate/"
DATASET_API_ENDPOINT = f"{BASE_URL}/datasets/api/dataset_generate/"

# Predictor list
PREDICTOR_LIST = [
    "arabidopsis_hvg20k",
    "corn_hvg10k",
    "rice_hvg20k",
    "soybean_hvg20k"
]

# Dataset list
DATASET_LIST = [
    "corn_SRP335_hvg10k",
    "rice_SRP286_hvg20k",
    "soybean_flowerbud_hvg20k",
    "SRP171_hvg20k",
    "SRP235_hvg20k",
    "SRP330_hvg20k",
    "tester"
]

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

for i in range(0, len(PREDICTOR_LIST)):
    payload = {
        "predictor_name": PREDICTOR_LIST[i],
        "predictor_filename": PREDICTOR_LIST[i],
        "predictor_public_flag": True
    }
    response = session.post(
        PREDICTOR_API_ENDPOINT,
        json=payload,
        headers=headers
    )
    if response.status_code == 200:
        print("API call successful!")
        print("Response:", response.json())
    else:
        print("API call failed!")
        print("Status Code:", response.status_code)
        print("Response:", response.text)

for i in range(0, len(DATASET_LIST)):
    payload = {
        "dataset_name": DATASET_LIST[i],
        "dataset_filename": DATASET_LIST[i],
        "dataset_public_flag": True
    }
    response = session.post(
        DATASET_API_ENDPOINT,
        json=payload,
        headers=headers
    )
    if response.status_code == 200:
        print("API call successful!")
        print("Response:", response.json())
    else:
        print("API call failed!")
        print("Status Code:", response.status_code)
        print("Response:", response.text)

# sudo cp -rf model_codebase/models/* core/uploads/predictors/chanye/
# sudo cp -rf model_codebase/public_data/* core/uploads/datasets/chanye/
# sudo cp -rf model_codebase/user_datasets/* core/uploads/datasets/chanye/

# ls -lha core/uploads/*/chanye/
# ls -lha model_codebase/public_data/ model_codebase/user_datasets/
