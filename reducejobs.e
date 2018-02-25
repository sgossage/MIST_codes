/n/home12/sgossage/.conda/envs/ANA_192/lib/python2.7/site-packages/matplotlib/__init__.py:878: UserWarning: text.fontsize is deprecated and replaced with font.size; please use the latter.
  warnings.warn(self.msg_depr % (key, alt_key))
Traceback (most recent call last):
  File "./reduce_jobs.py", line 47, in <module>
    reduce_jobs_utils.sort_histfiles(rawdirname)
  File "/n/home12/sgossage/MIST_codes/scripts/reduce_jobs_utils.py", line 196, in sort_histfiles
    mesa_hist_trim.trim_file(newhistfilename)
  File "/n/home12/sgossage/MIST_codes/scripts/mesa_hist_trim.py", line 40, in trim_file
    pagbind = np.where((logTeff > logTeff[logLHe[pagbind_tmp[0]]]+0.1) & ((starmass-ccoremass)/starmass[0] < 0.15))[0]
IndexError: only integers, slices (`:`), ellipsis (`...`), numpy.newaxis (`None`) and integer or boolean arrays are valid indices
