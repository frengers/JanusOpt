#!/usr/bin/env python

import sys, os, shutil, subprocess, time
import pandas as pd
import uuid
import jinja2 as jin
import argparse
import signal
import numpy as np


# Generate a input file
input_template = jin.Template('''
Input
Time = {{Time}}
dT = {{dT}}
tauc = {{tauc}}
taucWepp = {{taucWepp}}
lenzone = {{lenzone}}
n = {{n}}
nbare = {{nbare}}
Pmmphr = {{Pmmphr}}
tval = {{tval}}
''')

def get_args(argv):
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--pmmphr', dest='pmmphr')
    parser.add_argument('--tauc', dest='tauc')
    parser.add_argument('--taucwepp', dest='taucwepp')
    parser.add_argument('--maxtime', dest='maxtime')
    
    parser.set_defaults(pmmphr=40)
    parser.set_defaults(tauc=140)
    parser.set_defaults(taucwepp=6.2)
    parser.set_defaults(maxtime=900)
    
    return parser.parse_args(argv)

def create_params(args):
    
    params = {}
    params['dir_name'] = 'trial-'+str(uuid.uuid4())
    params['tval'] = 8000
    params['Pmmphr'] = args.pmmphr
    params['nbare'] = 0.022
    params['n'] = 0.05
    params['lenzone'] = 40
    params['taucWepp'] = args.taucwepp
    params['tauc'] = args.tauc
    params['dT'] = 1e-3
    params['Time'] = 10
    params['input_files'] = ['ForJanusOptimization.m', 'knickpointfunvDepSlope.m', 'profile.txt', 'xcell.mat', 'StrdxArray0.txt', 'StrElevArray0.txt']
    params['max_time'] = int(args.maxtime)
    params['wait_time'] = 15
    return params

def execute(params):
    
    directory = params['dir_name']
    cwd = os.getcwd()
    
    abs_path = os.path.abspath(directory)
    params['dir_name'] = abs_path
    
    # Create the directory
    if not os.path.exists(directory):
        os.mkdir(directory)
    
    # Create the input file
    tmp = input_template.render(params)
    with open(os.path.join(directory,'InputHack40.txt'),'w') as f:
        f.write(tmp)
        
    # Copy input files
    for file in params['input_files']:
        src = file
        dst = os.path.join(directory,file)
        shutil.copy(src, dst)
     
    # Run the model, maybe use subprocess
    os.chdir(directory)
    cmd = 'matlab -nodesktop -nosplash -nodisplay -r ForJanusOptimization > output'
    pro = subprocess.Popen(cmd, shell=True, preexec_fn=os.setsid)
    
    start = time.time()
    while pro.poll() is None:
        
        time.sleep(params['wait_time'])
        elapsed = (time.time() - start)
        #print "poll = ", pro.poll(), elapsed
        #if elapsed > 60:
            #print elapsed, pro.poll()
        if elapsed > params['max_time']:
            params['result'] = 1e9
            os.chdir(cwd)
            return params
            
    
    
    try:    
        os.killpg(pro.pid, signal.SIGTERM)
    except OSError, e:
        pass
    
    try:    
        # Need to compute objective function for search    
        last_z = np.genfromtxt('lastz.txt')
        profile = np.genfromtxt('profile.txt')
        assert len(last_z) == len(profile)
    
        #profile = np.linspace(min(last_z), max(last_z), len(last_z))
        tmp = (last_z - profile)
        res = np.dot(tmp,tmp)
        
        params['result'] = res
        
        with open('res.txt','w') as f:
            f.write(str(res)+'\n')
        
        with open('output_all.csv','w') as f:
            
            # HEader
            f.write('res,')
            f.write('tval,Pmmphr,n,lenzone,taucWepp,tauc,dT,Time,')
            last = len(last_z)-1
            for i,y in enumerate(last_z):
                if i == last:
                    f.write('y'+str(i)+'\n') 
                else:   
                    f.write('y'+str(i)+',')
            
            # Data
            f.write(str(res)+',')
            f.write(str(params['tval'])+',')
            f.write(str(params['Pmmphr'])+',')
            f.write(str(params['n'])+',')
            f.write(str(params['lenzone'])+',')
            f.write(str(params['taucWepp'])+',')
            f.write(str(params['tauc'])+',')
            f.write(str(params['dT'])+',')
            f.write(str(params['Time'])+',')
            for i,y in enumerate(last_z):
                if i == last:
                    f.write(str(y)+'\n')
                else:   
                    f.write(str(y)+',')
            
            with open('output.csv','w') as f:
            
                # HEader
                f.write('res,')
                f.write('tval,Pmmphr,n,lenzone,taucWepp,tauc,dT,Time\n')
                
                # Data
                f.write(str(res)+',')
                f.write(str(params['tval'])+',')
                f.write(str(params['Pmmphr'])+',')
                f.write(str(params['n'])+',')
                f.write(str(params['lenzone'])+',')
                f.write(str(params['taucWepp'])+',')
                f.write(str(params['tauc'])+',')
                f.write(str(params['dT'])+',')
                f.write(str(params['Time'])+'\n')
            
        os.chdir(cwd)        
        return params
    except:
        os.chdir(cwd)
        params['result'] = 1e9
        return params

def run(inputs):
    
  
    # Parse the command line args
    args = get_args(inputs)

    # Create params
    params = create_params(args)

    # Run
    result = execute(params)  
    
    # Cleanup
    return result

def collect_results(res):
    
    print 'Collecting results'
    
    frames = []
    for r in res:
        try:
            filename = os.path.join(r['dir_name'],'output.csv')
            if os.path.exists(filename):
                tmp = pd.read_csv(filename)
                frames.append(tmp)
            else:
                print filename
        except:
            pass
    
    data = pd.concat(frames)
    data.to_csv('final_output.csv', index=False)
    
    frames = []
    for r in res:
        try:
           filename = os.path.join(r['dir_name'],'output_all.csv')
           if os.path.exists(filename):
               tmp = pd.read_csv(filename)
               frames.append(tmp)
           else:
               print filename
        except:
            pass
    
    data = pd.concat(frames)
    data.to_csv('final_output_all.csv', index=False)

def main():        
        
    # Parse the command line args
    args = get_args(sys.argv[1:])

    # Create params
    params = create_params(args)

    # Run
    result = execute(params)  
    
    # Cleanup
    return result  
        
if __name__ == "__main__":
    
    res = main()
    print res




#os.system('ls -l trial-*')
