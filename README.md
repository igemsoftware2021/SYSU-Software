# Phoebus

[http://sysu-software.com](http://sysu-software.com)

## Environment

- OS: Ubuntu 18.04
- SDKs & Softwares:  Python 3.6.13, conda 4.10.1, Docker 20.10.9, Flask-Cors 3.0.10, Flask 1.1.2, nginx 1.16.1
- Important Notes: Please make sure port 4400 publicly available for viewing the 3d model of protein, port 5166 publicly available for accessing the Linker database, port 8097、8109 publicly available for accessing the File Server.

## Run Server


### Communication-Backend

#### How to run server?

You should install ```Python3``` first.

```bash
cd Communication-Backend
# Install requirements
pip3 install -r requirements.txt 

# Run Server
python3 community/manage.py runserver
# Then it will run on 127.0.0.1:8000

# When databases change(Not datas), please run below commands first
python3 community/manage.py makemigrations api
python3 community/manage.py migrate
# Then run the runserver command above
```
### MathModel-Backend

It's a small flask app to help with comunication between user and program.

#### How to run server?

You should install ```Python3``` first.

```bash
cd mathmodel-Backend
# Install requirements
pip install -r requirements.txt 

# Run Server
python main.py
# Then it will run on 127.0.0.1:8101
```
### Run Prediction Backend

The Prediction Backend runs for structure prediction algorithm, activity evaluator and sequence optimizer modules of our software.

```shell
cd Prediction-Backend
```

#### Prepare for Modified TrRosetta2

```shell
# create conda environment
conda env create -f Prediction-Backend.yml
conda activate Prediction-Backend

# download network weights [1.1G]
wget https://files.ipd.uw.edu/pub/trRosetta2/weights.tar.bz2
tar xf weights.tar.bz2

# download and install third-party software
./install_dependencies.sh

# download sequence and structure databases
# uniclust30 [46G]
wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz
mkdir -p UniRef30_2020_06
tar xf UniRef30_2020_06_hhsuite.tar.gz -C ./UniRef30_2020_06

# structure templates [8.3G]
wget https://files.ipd.uw.edu/pub/trRosetta2/pdb100_2020Mar11.tar.gz
tar xf pdb100_2020Mar11.tar.gz

# reinstall following packages for recent updates of these packages
conda install --force-reinstall biocore::blast-legacy=2.2.26
pip install git+https://github.com/songlab-cal/tape.git --upgrade --force-reinstall
```

Install PyRosetta [here](https://www.pyrosetta.org/home).

Verify the successful preparation through：

```shell
./run_pipeline.sh example/linker3_EAAAK_3.fasta example/linker3_EAAAK_3
```

Now you can fold proteins from scratch on your own machine！

#### Rebuild

```shell
g++ -O3 -o cadscore_1.1662/bin/voroprot2 src/*.cpp
```

#### Prepare for Backend

You should deploy redis first at port 6379.

```shell
sudo apt update
sudo apt install redis-server
sudo nano /etc/redis/redis.conf

. . .

# If you run Redis from upstart or systemd, Redis can interact with your
# supervision tree. Options:
#   supervised no      - no supervision interaction
#   supervised upstart - signal upstart by putting Redis into SIGSTOP mode
#   supervised systemd - signal systemd by writing READY=1 to $NOTIFY_SOCKET
#   supervised auto    - detect upstart or systemd method based on
#                        UPSTART_JOB or NOTIFY_SOCKET environment variables
# Note: these supervision methods only signal "process is ready."
#       They do not enable continuous liveness pings back to your supervisor.
supervised systemd
. . .

sudo systemctl restart redis.service
```

#### Get Started

```shell
cd prediction
celery worker -A app.tasks -l info
python run.py 
```

### Run Linker Database

The Linker database runs on the back-end, it is used to interact with the front-end and provides the linker data needed by the front-end.

#### Prepare for Linker Database

- First, you should deploy redis first at port 5166.
- Next, you should install apt or yum or binary packages as follows.

```shell
sudo apt install Flask==1.1.2
sudo apt install Flask-Cors==3.0.10
```

#### Get Started

```shell
cd linker/linker_database
python main.py
```

### Run File Server

#### Prepare for File Server

- First, you should Install the necessary dependencies.

```shell
sudo apt install build-essential
sudo apt install libtool
sudo apt install libpcre3 libpcre3-dev
sudo apt install zlib1g-dev
sudo apt install openssl
```

- Second, you should download the nginx binaries and the nginx upload module.

```shell
wget https://nginx.org/download/nginx-1.16.1.tar.gz
git clone https://github.com/vkholodkov/nginx-upload-module/
```

- Third, you should configure and install the nginx, the nginx installation path needs to be specified in the prefix argument.

```shell
cd nginx-1.16.1
./configure --with-compat --add-dynamic-module=../nginx-upload-module --prefix=[Nginx installation path]
make install
```

- then, you should modify the file server configuration file.

```
cd [Nginx installation path]/conf/nginx.conf
```

- You should add the following server configuration.

```shell
server {
        listen       8097 default_server;
        listen       [::]:8097 default_server;
        server_name  _;
        autoindex   on;

        location / {
            root   /home/SYSU/nginx/result_model;#Specifies the location where local files are stored by the server
            index  index.html index.htm;
        }

        error_page   500 502 503 504  /50x.html;
        location = /50x.html {
            root   html;
        }
    }

server {
    listen       8109 default_server;
    listen       [::]:8109 default_server;
    server_name  _;
    autoindex   on;

    location / {
        root   /home/SYSU/Opto;#Specifies the location where local files are stored by the server
        index  index.html index.htm;
    }

    error_page   500 502 503 504  /50x.html;
    location = /50x.html {
        root   html;
    }
}
```
#### Get Started

```shell
cd [Nginx installation path]/sbin
./nginx -c [Nginx installation path]/conf/nginx.conf
```
#### Access the file server for our project

You can access the file server already set up for our project by following the URL below. If you are interested in setting up a file server, you can download our data files at the following URL and use them to set up a file server.

[Linker PDB file server](http://20.106.156.143:8097/)

[Opto file server](http://20.106.156.143:8109/)

## Run Web Frontend App

The docker image `darkyzhou/sysu-software-fe` is bulit with the latest codes from the frontend application.

After installing Docker on the host, run this command:

```
docker run -d --restart unless-stopped --net host --name fe darkyzhou/sysu-software-fe:latest
```

After Docker finishes pulling images, the frontend application will be running. Visit `http://localhost` using web browsers to access the application.
