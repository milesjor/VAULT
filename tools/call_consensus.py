#!/usr/bin/env python3
# Chongwei 20200516
# bicwei@gmail.com

import os
import re
import subprocess
import glob
import shutil
import argparse
import multiprocessing
import sys


def get_argparse():
    parser = argparse.ArgumentParser(description='Just input the path to VAULT result folder '
                                                 '--> the folder contain (/snp /grouped_reads /umi analysis)')
    parser.add_argument('-s', '--save_path', type=str, required=True, help='path/to/VAULT_result_folder/')
    parser.add_argument('-t', '--thread', type=int, default='6', help='parallel worker [6]')
    parser.add_argument('-T', '--sub_thread', type=int, default='4', help='thread for every parallel worker [4]')
    parser.add_argument('--threshold', type=int, default='20', help='process UMI group with reads more than [20]')
    args = parser.parse_args()
    return args


# single worker, deprecated
def call_consensus():
    file_regex = './grouped_reads/perfect_umi/*/*.fastq'
    file_list = glob.glob(str(file_regex))
    thread = '50'
    save = './consensus'

    if os.path.isdir(save):
        print("folder %s exist, write file to the same folder" % save)
    else:
        os.mkdir(save)
        os.mkdir(save + '/medaka')
        os.mkdir(save + '/canu_draft')

    if len(file_list) >= 1:
        for file in file_list:
            file_name = re.split(r'/', file)[-1]
            file_name = re.split(r'.fastq', file_name)[0]
            read_count = re.split(r'_', file_name)[0]

            if int(read_count) >= 20:
                print("---> %s" % file)

                try:

                    cmd = """canu -d {save}/canu_draft/{file_name} -p {file_name} genomeSize=16500 -nanopore-raw {fastq} \
                                    useGrid=False maxThreads={thread} minInputCoverage=4 stopOnLowCoverage=4
                            """.format(file_name=file_name,
                                       fastq=file,
                                       thread=thread,
                                       save=save)
                    log = save + '/canu_draft/' + file_name + '.log'
                    with open(log, 'w') as log_file:
                        subprocess.run(cmd, shell=True, check=True, stdout=log_file, stderr=subprocess.STDOUT)

                except Exception as ex:
                    template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                    message = template.format(type(ex).__name__, ex.args)
                    print('\n' + message)
                    pass
                #
                # finally:
                #     break

                canu_draft = save + '/canu_draft/' + file_name + '/' + file_name + '.contigs.fasta'
                if os.path.isfile(canu_draft):
                    medaka_result = save + '/medaka/' + file_name

                    try:
                        cmd = """medaka_consensus -i {fastq} -d {canu_draft} -o {medaka_result} -t {thread} -m r941_min_high_g344
                                # &> {medaka_result}/{file_name}.log
                                """.format(fastq=file,
                                           canu_draft=canu_draft,
                                           medaka_result=medaka_result,
                                           thread=thread,
                                           file_name=file_name)
                        log = save + '/medaka/' + file_name + '.log'
                        with open(log, 'w') as log_file:
                            subprocess.run(cmd, shell=True, check=True, stdout=log_file, stderr=subprocess.STDOUT)

                    except Exception as ex:
                        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
                        message = template.format(type(ex).__name__, ex.args)
                        print('\n' + message)
                        pass

                    medaka_consensus = medaka_result + '/consensus.fasta'

                    if os.path.isfile(medaka_consensus):
                        new_name = save + '/' + file_name + '.consensus.fasta'
                        shutil.move(medaka_consensus, new_name)
                    else:
                        print("<--- medaka failed for %s" % file)

                else:
                    print("<--- canu failed for %s" % file)

    else:
        print("No fastq found in %s" % file_regex)


# multiprocessing
def single_cmd(file, file_name, thread, save):
    # print("Proceeding %s" % file)

    try:

        cmd = """canu -d {save}/canu_draft/{file_name} -p {file_name} genomeSize=16500 -nanopore-raw {fastq} \
                 useGrid=False maxThreads={thread} minInputCoverage=4 stopOnLowCoverage=4""".format(file_name=file_name,
                                                                                                    fastq=file,
                                                                                                    thread=thread,
                                                                                                    save=save)
        log = save + '/canu_draft/' + file_name + '.log'
        with open(log, 'w') as log_file:
            subprocess.run(cmd, shell=True, check=True, stdout=log_file, stderr=subprocess.STDOUT)

    except Exception as ex:
        template = "WARNING!  An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)
        pass
    #
    # finally:
    #     break

    canu_draft = save + '/canu_draft/' + file_name + '/' + file_name + '.contigs.fasta'
    if os.path.isfile(canu_draft):
        medaka_result = save + '/medaka/' + file_name

        try:
            cmd = """medaka_consensus -i {fastq} -d {canu_draft} -o {medaka_result} -t {thread} -m r941_min_high_g344
                    """.format(fastq=file,
                               canu_draft=canu_draft,
                               medaka_result=medaka_result,
                               thread=thread,
                               file_name=file_name)
            log = save + '/medaka/' + file_name + '.log'
            with open(log, 'w') as log_file:
                subprocess.run(cmd, shell=True, check=True, stdout=log_file, stderr=subprocess.STDOUT)

        except Exception as ex:
            template = "WARNING!  An exception of type {0} occurred. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print('\n' + message)
            pass

        medaka_consensus = medaka_result + '/consensus.fasta'

        if os.path.isfile(medaka_consensus):
            new_name = save + '/' + file_name + '.consensus.fasta'
            shutil.move(medaka_consensus, new_name)
        else:
            print("WARNING!  medaka failed for %s\n" % file)

    else:
        print("WARNING!  canu failed for %s\n" % file)


def parallel_run(args):
    file_regex = args.save_path + '/grouped_reads/perfect_umi/*/*.fastq'
    file_list = glob.glob(str(file_regex))
    thread = args.thread
    sub_thread = args.sub_thread
    save = args.save_path + '/consensus'

    if os.path.isdir(save):
        print("WARNING!  folder %s exist, write file to the same folder" % save)
    else:
        os.mkdir(save)
        os.mkdir(save + '/medaka')
        os.mkdir(save + '/canu_draft')

    if len(file_list) >= 1:
        p = multiprocessing.Pool(thread)
        for file in file_list:
            file_name = re.split(r'/', file)[-1]
            file_name = re.split(r'.fastq', file_name)[0]
            read_count = re.split(r'_', file_name)[0]

            if int(read_count) >= args.threshold:
                p.apply_async(single_cmd, args=(file, file_name, sub_thread, save))

        p.close()
        p.join()

    else:
        sys.stderr.write('\nERROR!  No fastq file detected in %s/grouped_reads/perfect_umi/*/*.fastq\n'
                         'ERROR!  Please check input folder\n' % args.save_path)
        sys.exit(1)


def check_tool(tool):
    """Check whether `tool` is on PATH and marked as executable."""
    return shutil.which(tool) is not None


def consensus_main(args):
    if check_tool("canu") is not True:
        sys.exit("ERROR! Executable canu is not found!\n"
                 "ERROR! Please install canu -> `conda install canu`")
    if check_tool("medaka") is not True:
        sys.exit("ERROR! Executable medaka is not found!\n"
                 "ERROR! Please install medaka -> `conda install medaka`")

    print("\n-------------- received arguments -------------")
    print("VAULT folder      :   %s" % args.save_path)
    print("parallel worker   :   %s" % args.thread)
    print("thread per worker :   %s" % args.sub_thread)
    print("read number threshold :   %s" % args.threshold)
    print("---> Start analyzing ...\n")
    parallel_run(args)
    print("------------ Jobs done! ------------\n")
    sys.exit(0)


if __name__ == '__main__':
    args = get_argparse()
    consensus_main(args)
