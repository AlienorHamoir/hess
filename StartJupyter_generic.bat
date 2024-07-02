set username=alienor
set container=cgehbauer/jupyter_radiance_eplus:v3
start "" http://127.0.0.1:8889/
docker run -p 127.0.0.1:8889:8888 -v data:/home/%username%/data -v C:\Users\%username%:/home/%username% -it --rm %container% bash -c "cd /home/%username% && jupyter notebook --ip=0.0.0.0 --allow-root --no-browser --NotebookApp.token=''"
PAUSE