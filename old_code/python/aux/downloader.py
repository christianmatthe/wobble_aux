import os
import shutil
import requests

def download_folder(url, folder_name):
    local_filename = url.split('/')[-1]
    path = os.path.join("/{}/{}".format(folder_name, local_filename))
    print(path)
    user, password = 'cmatthe', 'rectangular box apparatus'
    with requests.get(url, stream=True, auth=(user, password)) as r:
        with open(path, 'wb') as f:
            shutil.copyfileobj(r.raw, f)

    return local_filename

if __name__ == "__main__":
    site_url= "https://svn.mpia.de/trac/gulli/carmenes-science/login"
    url = "https://svn.mpia.de/trac/gulli/carmenes-science/browser/data/trunk/serval/CARM_VIS/J11302%2B076/J11302%2B076.rvc.dat"
    #url = "http://carmenes.cab.inta-csic.es/gto/getDataAction.action?id=car-20171226T19h27m20s-sci-gtoc-vis_A.fits"
    
    user, password = 'cmatthe', 'rectangular box apparatus'
    #resp = requests.get(url, auth=(user, password))
    
    #path = "/data/cmatthe/test4.dat"
    #with requests.get(url, stream=True, auth=(user, password)) as r:
        #with open(path, 'wb') as f:
            #shutil.copyfileobj(r.raw, f)
   

    #download_folder("https://svn.mpia.de/trac/gulli/carmenes-science/browser/data/trunk/serval/CARM_VIS/", "J22020-194") # NOTE needs work and apparently more than my current coding prowess

    #### wget attempt
    '''
    # Log in to the server.  This only needs to be done once.
    wget --save-cookies cookies.txt \
    --keep-session-cookies \
    --post-data 'user=foo&password=bar' \
    --delete-after \
    http://server.com/auth.php

    # Now grab the page or pages we care about.
    wget --load-cookies cookies.txt \
    http://server.com/interesting/article.php
    
    
    #os.system("wget --user ={0} --ask-password {1}".format(user, url))
    post_string = 'user={0}&password={1}'.format()
    os.system(wget --save-cookies cookies.txt \
    --keep-session-cookies \
    --post-data 'user={0}&password={1}'.format \
    --delete-after \
    site_url)
    os.system("wget --user ={0} --password = {1} {2}".format(user, password, url))
    '''
    
    #persistent cookies with requests
    #path = "/data/cmatthe/test4.dat"
    #s = requests.Session()
    #s.get(site_url)
    #s.post(site_url, data={'user': user, 'password': password})
    ##s.get(url)
    #with s.get(url, stream=True) as r:
        #with open(path, 'wb') as f:
            #shutil.copyfileobj(r.raw, f)
            
            
    #from subprocess import Popen, PIPE
    #from time import sleep
    #from fcntl import fcntl, F_GETFL, F_SETFL
    #from os import O_NONBLOCK, read

    ## run the shell as a subprocess:
    #p = Popen(['python', 'shell.py'],
            #stdin = PIPE, stdout = PIPE, stderr = PIPE, shell = False)
    ## set the O_NONBLOCK flag of p.stdout file descriptor:
    #flags = fcntl(p.stdout, F_GETFL) # get current p.stdout flags
    #fcntl(p.stdout, F_SETFL, flags | O_NONBLOCK)
    ## issue command:
    #p.stdin.write('command\n')
    ## let the shell output the result:
    #sleep(0.1)
    ## get the output
    #while True:
        #try:
            #print read(p.stdout.fileno(), 1024),
        #except OSError:
            ## the os throws an exception if there is no data
            #print '[No more data]'
            #break
        
        
        
    ###
    import subprocess
    #command = 'echo a; echo b'
    
    command = "wget --save-cookies cookies.txt \
    --keep-session-cookies \
    --post-data 'user={0}&password={1}' \
    --delete-after \
    {2} ;wget --load-cookies cookies.txt \
    {3}".format(user, password, site_url, url)
    
    
    process = subprocess.Popen(command,stdout=subprocess.PIPE, shell=True)
    proc_stdout = process.communicate()[0].strip()
    print(proc_stdout)
    
    
