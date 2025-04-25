## Visualizing Anvi'o Products 

To do this, you will need two terminals open simultaneously. In the first terminal window you will set up Anvi'o. In the second, you will create an SSH tunnel.


### Logging in to the teaching cluster 
On the first terminal, you will have to sign into the teaching cluster and request an interactive node. 

```
ssh userid@teach.gacrc.uga.edu
interact mem=12G
```

Then you will have to obtain the IP address of the host server.

```
hostname -i 
```


### Setting up SSH Tunnel
Now, on the second terminal, you will have to set up an SSH Tunnel. You will do this by specifying two flags `-N` `-L` to redirect the output to the port on your local computer. A port is a connection that can handle incoming and outgoing request. Think of it as specifying the road that the program or computer should use. 

```
ssh -N -L PORT:HOST:PORT user@xfer.gacrc.uga.edu
```

Just replace PORT with a PORT number (a common one is 8080) and change the HOST to the IP Address you saved earlier. If you get a message the the line is busy, change the port number. 

for example, if the IP address is 88.88.88.888 then you would run the following code. 

```
ssh -N -L 8080:88.88.88.888:8080 userid@teach.gacrc.uga.edu
```

Note that when you create the ssh tunnel it will look as though the screen is frozen. DO NOT PANIC - that means that the connection is working. 

### Setting up Anvi'o
Now, let's go back to the first terminal and load our conda environment.

```
module load Miniconda3
source activate /home/ad14556/conda-env/anvio/
```

Now, let's navigate to the location of our anvi'o outputs.

```
cd /work/mars8180/instructor_data/metagenomic-datasets/nematode-microbiome/results/18-anvio-short-reads/epacanthion.1
```

Now, you can run the `anvi-interactive` command to visualize your data. You will have to use -P to specify the port that Anvi'o should use. Following the example above, if you are using port 8080, you will run the command like this: 

```
anvi-interactive -p profile/PROFILE.db -c epacanthion.1.db --server-only -P 8080
```

Now, start a browser on your computer (preferably Chrome) and type the following address `http://localhost:8080`

