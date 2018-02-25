/n/home12/sgossage/.conda/envs/ANA_192/lib/python2.7/site-packages/matplotlib/__init__.py:878: UserWarning: text.fontsize is deprecated and replaced with font.size; please use the latter.
  warnings.warn(self.msg_depr % (key, alt_key))
Traceback (most recent call last):
  File "./reduce_jobs.py", line 71, in <module>
    make_eeps_isos.make_eeps_isos(runname, basic=False)
  File "/n/home12/sgossage/MIST_codes/scripts/make_eeps_isos.py", line 52, in make_eeps_isos
    make_iso_input_file.make_iso_input_file(runname, "eeps", basic, custom_path=custom_path)
  File "/n/home12/sgossage/MIST_codes/scripts/make_iso_input_file.py", line 85, in make_iso_input_file
    fehstr, feh, afestr, afe, vvcrit, extra = dirname.split('_')
ValueError: too many values to unpack
