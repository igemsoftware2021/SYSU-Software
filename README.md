# Phoebus

\<domain name\>

## Environment

- OS: Ubuntu 18.04
- SDKs & Softwares:  Python 3.6.13, conda 4.10.1, Docker 20.10.9
- Important Notes: Please make sure port 4400 publicly available for viewing the 3d model of protein.

## Run Server

### Run Prediction Backend

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



### Run Web Frontend App

The docker image `darkyzhou/sysu-software-fe` is bulit with the latest codes from the frontend application.

After installing Docker on the host, run this command:

```
docker run -d --restart unless-stopped --net host --name fe darkyzhou/sysu-software-fe:latest
```

After Docker finishes pulling images, the frontend application will be running. Visit `http://localhost` using web browsers to access the application.
