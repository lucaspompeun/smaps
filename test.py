def bowtie2(read1, reference, project, read2 = None, N='1', L='22', threads='16'):
    print("bowtie2...")
    """
    Run bowtiw2 and returns sam file
    """
    out = project + 'bowtie/'
    database = out + "database"
    if not os.path.exists(out):
        os.mkdir(out)
    os.system("bowtie2-build " + reference + " " + database + " >  " + out + 'database.log')
    
    cmd = "bowtie2 -p " + threads + " -x " + database 
    if read2:
        cmd += " -1 " + read1 + " -2 " + read2

    else:
        cmd += " -U " + read1
        
    cmd += " -N " + N + " -L " + L + " -S " + out + "output.sam > " + out + 'execution.log'

    write_file(out + 'comandline.txt', cmd)

    os.system(cmd)

    return out + 'output.sam'