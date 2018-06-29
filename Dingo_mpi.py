#!/usr/bin/env python

"""

Project_Name: main, File_name: master_mpi
Aufthor: kalabharath, Email: kalabharath@gmail.com
Date: 13/04/15 , Time:10:05 AM

One master mpi and one master_search to eliminate redundancy across four files
should be called as a part of the sequence of smotif assembly files

"""

# sys.path.append('../../main/')
import argparse
import time
import traceback
from   mpi4py import MPI

import ranking.SmotifRanking as srank
from ranking.NoeStageRank import *
import smotif_search as msearch
import utility.masterutil as mutil
import utility.stage2_util as util
import utility.io_util as io

# Define MPI message tags

tags = mutil.enum('READY', 'DONE', 'EXIT', 'START')
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
name = MPI.Get_processor_name()
status = MPI.Status()


def killall(processes):
    """
    Kill all the subprocess when requested
    :param processes:
    :return: True or False
    """
    count = 0
    while True:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            comm.send(None, dest=source, tag=tags.EXIT)
            count += 1
        if count == processes -1:
            break
    return True


# ************************************  Define cmd line argument parser ************************************

parser = argparse.ArgumentParser(description='DINGO-PCS Master MPI process that manages all jobs.')
parser.add_argument('--stage', type=int, help='specify the stage of  the Smotif assembly')
parser.add_argument('--numhits', type=int, help='Top number of hits to be selected from previous assembly')
parser.add_argument('--infile', type=int, help='specify the top_hits file')
args = parser.parse_args()
# ************************************ Define cmd line argument parser ************************************

# Rank '0' specifies the master process

if rank == 0:

    # ************************************ Extract top hits from previous stage ************************************

    if args.stage == 1:
        tasks = mutil.getRunSeq()  # there are no hits to extract if it is the 1st stage

    else:
        # for 2nd stage and beyond
        try:
            try:
                # Restart from the already assembled top hits if possible
                tasks, sse_index = util.start_top_hits(args.numhits, args.stage, args.infile)
                print "here", tasks, sse_index
            except:
                # Assemble top hits from the previously generated hits
                print "XXX:", args.numhits, args.stage
                tasks, sse_index = srank.getRunSeq(args.numhits, args.stage)  # TODO change this for new alt_smotifs
                # tasks, sse_index = srank.getRunSeqAlt(args.numhits, args.stage, args.infile)
        except:
            # print what went wrong and terminate the slave processes
            traceback.print_exc()
            print "Couldn't extract top hits within the specified cutoffs: Exiting..."
            killall(size)
            exit()

    # ************************************ Generate and distribute job index array ************************************

    stime = time.time()

    try:
        if len(tasks):
            pass
    except:
        print "killing all processes!"
        killall(size)
        exit()

    # print tasks, len(tasks) # this will be the new tasks
    task_index = 0  # control the number of processes with this index number
    finished_task = 0
    num_workers = size - 1  # 1 processor is reserved for master.
    closed_workers = 0  # control the workers with no more work that can be assigned

    print ("Master starting with {} workers".format(num_workers))
    total_data = []

    while closed_workers < num_workers:
        # Manage/distribute all processes in this while loop
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # worker process is ready, send some task to do.
            if task_index < len(tasks):
                comm.send([tasks[task_index], args.stage, args.infile], dest=source, tag=tags.START)
                task_index += 1  # increment its
            else:
                # everything is done, send exit signal
                comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.DONE:
            # take the result from the worker
            if data:
                for hit in data:
                    total_data.append(hit)
            ctime = time.time()
            elapsed = ctime - stime
            finished_task += 1
            print "Finishing..", finished_task, "of", len(tasks), "Smotifs, Elapsed", round((elapsed) / (60), 2), "mins"
        elif tag == tags.EXIT:
            closed_workers += 1

    # consolidate top_hits and dump files here
    if args.stage == 1:
        tasks, sse_index = srank.getRunSeq(args.numhits, args.stage)  # TODO change this for new alt_smotifs
        exit()
    print "Total number of hits  found are : ",len(total_data)
    # ranked_data = rank_assembly(total_data, args.numhits)
    ranked_data = rank_assembly_with_clustering(total_data, args.numhits)
    print len(ranked_data)
    if args.stage == 1:
        sse_index = 0
    io.dumpGzipPickle(str(sse_index) + "_tophits.gzip", ranked_data)
    # Rename temprary files
    util.rename_pickle(sse_index)
    print "All Done, Master exiting"
    exit()


# On the worker processes
else:

    while True:  # initiate infinite loop
        comm.send(None, dest=0, tag=tags.READY)
        # Signal the master process that you are READY

        task = comm.recv(source=0, tag=MPI.ANY_SOURCE, status=status)
        tag = status.Get_tag()
        if tag == tags.START:
            result = False
            if args.stage == 1:
                result = msearch.S1SmotifSearch(task)
            else:
                result = msearch.sXSmotifSearch(task)

            comm.send(result, dest=0, tag=tags.DONE)
        elif tag == tags.EXIT:
            # break the infinite loop because there is no more work that can be assigned
            break

    # Signal EXIT to the master process
    comm.send(None, dest=0, tag=tags.EXIT)
