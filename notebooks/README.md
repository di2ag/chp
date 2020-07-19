instructions for port-forwarding to allow use of jupyter notebooks. Make sure you are on a compute node!!!

1. create a tmux session with:
 
tmux new -s <session_name>

1.5 make sure you add your virtual environment to jupyter notebook via ipykernel. Run the following to do so:

python -m ipykernel install --user --name=<env_name>

2. launch a headless jupyter notebook session with:

jupyter notebook --no-browser --port=8888

3. jupyter notebook will print out information in the terminal. Collect the token id. This should look something like,
"http://localhost:8888/?token=f7c6c110319c71af763b3fc52ba4e0a5f57d15f9c100e6b1", where the token you need to copy is,
"f7c6c110319c71af763b3fc52ba4e0a5f57d15f9c100e6b1" 

4. launch some ssh client software from your local machine (windows cmd is fine) and port forward using

ssh -L 8888:localhost:8888 <user_name>@intuition.thayer.dartmouth.edu ssh -L 8888:localhost:8888 -N <compute_node_id>
where <user_name> is your login id for intution and <compute_node_id> is the name of the compute node (i.e. c-dell-m630-0-16)

5. Open a browser on your local machine and paste:

localhost:8888

6. You will be asked to enter a token. Paste the token you copied from step 3.

7. When you open a notebook you will have to select your appropriate virtual environment. From the menu select:
kernal
    -> Change kernel
        -> <env_name>

You're all set to go!!! 

#####################################
#            Important              #
#####################################
When you leave jupyter notebook make sure to kill your ssh tunnel. This should be done automatically,
but occasionally is not. To check if it was killed, mode to the intuition head node (NOT the compute node)
and run

lsof -ti:8888

if nothing was returned, the tunnel was succesfully closed, otherwise if there was terminal output run

lsof -ti:8888 | xargs kill -9
