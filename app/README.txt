# Only needed for headless Xserver
sudo apt update
sudo apt install xvfb xauth xfonts-base
Xvfb :99 &
export DISPLAY=:99

pip3 install -r requirements.txt
python3 -m streamlit run ./app.py
