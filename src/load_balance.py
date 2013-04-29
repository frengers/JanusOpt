#!/usr/bin/env python
from mpi4py import MPI
import math
import os, sys
import subprocess
from time import sleep
import socket
#import tau

# Create some enumerated types
def enum(**enums):
    return type('Enum', (), enums)

tags = enum(WORK=1, WORK_DONE=2, KILL=3)


def simulator(cmd):
    
    #print cmd
    #cmd = 'hostname; echo 10; pwd'
    #print cmd
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    host = socket.gethostname()
   # print host, rank, "start"
    #os.system('hostname')
    #print cmd
    os.system(cmd)
    #print host, rank, "end"
    # try:
  #       retcode = subprocess.Popen(cmd)
  #       #retcode = call(cmd, shell=True)
  #       
  #       retcode.wait()
  #       
  #       if not retcode.poll():
  #           print >>sys.stderr, "Child was TERMINATED by signal", retcode.poll(), host, rank
  #       else:
  #           print >>sys.stderr, "SUCCESS", retcode.poll(), host, rank
  #   except OSError as e:
  #       print >>sys.stderr, "FAILED", e, host, rank
        
   #print os.uname()[1], str(cmd)
    #output = subprocess.check_output([cmd] , shell=True)
    return True

def create_queue(filename):
    import os
    file_in = open(os.path.join(os.getcwd(), filename),'r') 
    line = file_in.readline().rstrip('\n')
    work_queue = []
    while line:
        #print line
        work_queue.append(line)
        line = file_in.readline().rstrip('\n')
    file_in.close()
    return work_queue
    
def create_queue_args(args_jobs, args_time):
   
    work_queue = []
    for i in range(args_jobs):
        work_queue.append(args_time)
        
    return work_queue

#not used right now...
def kill_switch(comm, rank, size):
    if rank == 0:
        requests = [MPI.REQUEST_NULL] * (size-1)
        for r in range(1,size):
            index = r-1
            data = 0
            requests[index] = comm.isend(data, r, tags.KILL)
        MPI.Request.Waitall(requests)        

def load_balance(work_queue, groups, comm, size):
    
    print groups
    #print groups
    requests = [MPI.REQUEST_NULL] * (size-1)
   
    start_index = int(0)
    end_index = int(start_index + groups[0])
    
    for r in range(0,size-1):
        end_index = int(start_index + groups[r])
        #print work_queue[start_index:end_index]
        tmp = work_queue[start_index:end_index]
        requests[r] = comm.isend(tmp, r+1, tags.WORK) 
        start_index = end_index
    
    MPI.Request.Waitall(requests)    
    
    # Master does this work
    end_index = int(start_index + groups[-1])        
    for i in work_queue[start_index:end_index]:
        simulator(i)
    
   
       
    # Kill switch... wait for each rank to finish...
    for r in range(1,size):
        data = 0
        comm.send(data, r, tags.KILL)

                
def get_args(argv, size):
    import argparse 
    parser = argparse.ArgumentParser(
                    description='A python load balancer',
                    epilog='monte.lunacek@colorado.edu')
                    
    parser.add_argument('-f', '--file', help='Command lines to parse')
    
    return parser.parse_args(argv)
    
def worker(comm, rank, size):
   
    while True:
        status = MPI.Status()
        msg = comm.Probe(0, MPI.ANY_TAG, status=status)
        if status.Get_tag() == tags.WORK:
            data = comm.recv(source=0, tag=tags.WORK)
            for i in data:
                simulator(i)
        
        if status.Get_tag() == tags.KILL:
            data = comm.recv(source=0, tag=tags.KILL)
            done = True
            break   

def split_jobs(size, N):
    
    num_groups = int(size)
    small_group = math.floor(N/float(num_groups))
    large_group = math.ceil(N/float(num_groups))
    if small_group == large_group:
    	num_small_group = N/small_group
    	num_large_group = 0
    else:
    	num_large_group = (N - num_groups*small_group)/(large_group-small_group)
    	num_small_group = num_groups - num_large_group

    #print N, num_large_group, large_group, num_small_group, small_group
    group_array = []
    for i in range(int(num_large_group)):
         group_array.append(large_group)
         
    for i in range(int(num_small_group)):
         group_array.append(small_group)

    return group_array

def get_work_queue(comm, rank, size):
    
    import sys
    args = get_args(sys.argv[1:], size)
    if not args.file:
        return None
        
    work_queue = create_queue(args.file)
    return work_queue

def master(comm, rank, size):
    
    import sys
    args = get_args(sys.argv[1:], size)
    
    if not args.file:
        # Kill ranks
        for r in range(1,size):
            data = 0
            comm.send(data, r, tags.KILL)
        return False
    
    work_queue = get_work_queue(comm, rank, size)
    num_jobs = len(work_queue)
    groups = split_jobs(size, num_jobs)
   
    load_balance(work_queue, groups, comm, size)
    
    return True

def main():
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    # num_jobs = None
 #    if rank == 0:
 #        work_queue = get_work_queue(comm, rank, size)
 #        num_jobs = len(work_queue)
 #    
 #    num_jobs = comm.bcast(num_jobs, root=0)
 #    
 #    groups = split_jobs(size, num_jobs)
 #    print groups
    
    # Send each rank it's jobs

    if rank == 0:
        master(comm, rank, size)
    else:
        worker(comm, rank, size)



if __name__ == '__main__':
    main()

#tau.run('main()')   
    
    
    
    
    
    
    
    
    
    
    
    
