import urllib.request
import os

# List generated with AI
# the files from https://github.com/rub-hgi/shorter_linear_slps_for_mds_matrices/tree/master/matrices_bp_format
txt_files = [
    "ACISP_SarSye17_8x8_4.txt",
    "ACISP_SarSye17_8x8_8.txt",
    "AES.txt",
    "Anubis.txt",
    "C_BeiKraLea16_4x4_4.txt",
    "C_BeiKraLea16_4x4_8.txt",
    "C_BeiKraLea16_8x8_8.txt",
    "Clefia_M0.txt",
    "Clefia_M1.txt",
    "FSE_LiWang16_4x4_4.txt",
    "FSE_LiWang16_4x4_4_2.txt",
    "FSE_LiWang16_4x4_8.txt",
    "FSE_LiWang16_4x4_8_2.txt",
    "FSE_LiWang16_i_4x4_4.txt",
    "FSE_LiWang16_i_4x4_8.txt",
    "FSE_LiWang16_i_4x4_8_2.txt",
    "FSE_LiuSim16_4x4_4.txt",
    "FSE_LiuSim16_4x4_8.txt",
    "FSE_LiuSim16_8x8_8.txt",
    "FSE_SKOP15_4x4_4.txt",
    "FSE_SKOP15_4x4_8.txt",
    "FSE_SKOP15_8x8_4.txt",
    "FSE_SKOP15_8x8_8.txt",
    "FSE_SKOP15_i_4x4_4.txt",
    "FSE_SKOP15_i_4x4_8.txt",
    "FSE_SKOP15_i_8x8_4.txt",
    "FSE_SKOP15_i_8x8_8.txt",
    "Fox_Mu4.txt",
    "Fox_Mu8.txt",
    "Grostl.txt",
    "Joltik.txt",
    "Khazad.txt",
    "MIDORI.txt",
    "M_4_4.txt",
    "M_4_8.txt",
    "M_8_4.txt",
    "M_8_8.txt",
    "M_i_4_8.txt",
    "M_i_8_4.txt",
    "M_i_8_8.txt",
    "PRIDE_L_0.txt",
    "PRIDE_L_1.txt",
    "PRIDE_L_2.txt",
    "PRIDE_L_3.txt",
    "PRINCE_M_0.txt",
    "PRINCE_M_1.txt",
    "QARMA128.txt",
    "QARMA64.txt",
    "SKINNY.txt",
    "SmallScale_AES.txt",
    "ToSC_SarSye16_4x4_4.txt",
    "ToSC_SarSye16_4x4_8.txt",
    "ToSC_SarSye16_i_4x4_4.txt",
    "ToSC_SarSye16_i_4x4_8.txt",
    "Twofish.txt",
    "Whirlpool.txt",
    "Whirlwind_M0.txt",
    "Whirlwind_M1.txt",
    "ePrint_JeaPeySim_4x4_4.txt",
    "ePrint_JeaPeySim_4x4_8.txt",
    "ePrint_JeaPeySim_i_4x4_4.txt",
    "ePrint_JeaPeySim_i_4x4_8.txt",
    "ePrint_JeaPeySim_i_8x8_8.txt",
]

os.makedirs("benchmarks/crypt", exist_ok=True)
output_dir = "benchmarks/crypt/dataset"
os.makedirs(output_dir, exist_ok=True)

url_base = "https://raw.githubusercontent.com/rub-hgi/shorter_linear_slps_for_mds_matrices/refs/heads/master/matrices_bp_format/"

for file in txt_files:
    download_url = url_base + file
    destination = os.path.join(output_dir, file)
    print(f"Downloading {file}...")
    urllib.request.urlretrieve(download_url, destination)

print("All files downloaded successfully!")
