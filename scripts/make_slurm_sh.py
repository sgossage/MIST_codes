"""

Generates the Odyssey cluster SLURM file for each model.

Args:
    inlistname: name of the inlist
    inlistdir: name of the inli
    runbasefile: the name of the template shell script
    
Returns:
    the name of the SLURM file

"""

def make_slurm_sh(inlistname, inlistdir, runbasefile, partition):
    
    #Read the contents of the base file 
    infile = open(runbasefile, 'r')
    infile_contents = infile.read()
    infile.close()

    #Replace the keys with appropriate values
    runname = inlistname.strip(".inlist")
    replaced_contents = infile_contents.replace('<<RUNNAME>>', runname)
    replaced_contents = replaced_contents.replace('<<DIRNAME>>', inlistdir)
    replaced_contents = replaced_contents.replace('<<PARTITION>>', partition)

    # Determine number of cores to use; only two partitions are defined -- system specific.
    if partition == 'conroy':
        num_cores = '8'
    elif partition == 'conroy-intel':
        num_cores = '4'
    else:
        print("Partition \"{:s}\" has an unspecified num_cores in make_slurm_sh.py.")
        return ""

    replaced_contents = replaced_contents.replace('<<CORES>>', num_cores)

    #Find the mass of the model to assign appropriate runtime
    massval = int(runname.split('M')[0])/100.0
    if (massval <= 0.6):
        # Originally 30:00:00, added 20 hrs
        replaced_contents = replaced_contents.replace('<<RUNTIME>>', '50:00:00')
    else:
        # Originally 120:00:00, added 20 hrs
        replaced_contents = replaced_contents.replace('<<RUNTIME>>', '170:00:00')
    runfile = runname+'_run.sh'

    #Write the new shell script
    outfile = open(runfile, 'w+')
    outfile.write(replaced_contents)
    outfile.close()

    return runfile
