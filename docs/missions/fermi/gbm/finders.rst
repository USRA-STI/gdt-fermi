.. _gbm-finders:
.. |TriggerFtp| replace:: :class:`~gdt.missions.fermi.gbm.finders.TriggerFtp`
.. |ContinuousFtp| replace:: :class:`~gdt.missions.fermi.gbm.finders.ContinuousFtp`

**********************
Fermi GBM Data Finders
**********************
A natural question may be: "Where do I find the data I need?" Well, you're in 
luck, because this will show you how to find the data you seek. GBM Data is 
hosted publicly on the HEASARC FTP server via the Fermi Science Support Center, 
and the data are stored in a consistent directory structure. But instead of 
having to navigate a winding maze of FTP directories, we provide a couple of 
classes built to retrieve the data you want. First, you need to decide if you 
want trigger data (say, from a GRB) or continuous data. 

Finding Triggered GBM Data
==========================

Let's start with trigger data, and assume you know the trigger number you're 
interested in (190114873):

    >>> from gdt.missions.fermi.gbm.finders import TriggerFtp
    >>> trig_finder = TriggerFtp('190114873')
    <TriggerFtp: 190114873>
    >>> trig_finder.num_files
    122

We don't really care about the directory structure, we just want the data. So 
this quickly gets us to the directory we need. There are 122 files associated 
with this trigger. Say we want CSPEC data. is there CSPEC available?

    >>> trig_finder.ls_cspec()
    ['glg_cspec_b0_bn190114873_v00.pha',
     'glg_cspec_b1_bn190114873_v00.pha',
     'glg_cspec_n0_bn190114873_v00.pha',
     'glg_cspec_n1_bn190114873_v00.pha',
     'glg_cspec_n2_bn190114873_v00.pha',
     'glg_cspec_n3_bn190114873_v00.pha',
     'glg_cspec_n4_bn190114873_v00.pha',
     'glg_cspec_n5_bn190114873_v00.pha',
     'glg_cspec_n6_bn190114873_v00.pha',
     'glg_cspec_n7_bn190114873_v00.pha',
     'glg_cspec_n8_bn190114873_v00.pha',
     'glg_cspec_n9_bn190114873_v00.pha',
     'glg_cspec_na_bn190114873_v00.pha',
     'glg_cspec_nb_bn190114873_v00.pha']

Great! There's a full complement of CSPEC data. How about responses for the 
CSPEC data?

    >>> trig_finder.ls_rsp(cspec=True, ctime=False)
    ['glg_cspec_b0_bn190114873_v02.rsp',
     'glg_cspec_b1_bn190114873_v02.rsp',
     'glg_cspec_n0_bn190114873_v02.rsp',
     'glg_cspec_n1_bn190114873_v02.rsp',
     'glg_cspec_n2_bn190114873_v02.rsp',
     'glg_cspec_n3_bn190114873_v02.rsp',
     'glg_cspec_n4_bn190114873_v02.rsp',
     'glg_cspec_n5_bn190114873_v02.rsp',
     'glg_cspec_n6_bn190114873_v02.rsp',
     'glg_cspec_n7_bn190114873_v02.rsp',
     'glg_cspec_n8_bn190114873_v02.rsp',
     'glg_cspec_n9_bn190114873_v02.rsp',
     'glg_cspec_na_bn190114873_v02.rsp',
     'glg_cspec_nb_bn190114873_v02.rsp']

What if we want to move on to another trigger? You don't have to create a new 
|TriggerFTP| object, you can just used ``cd()``:

    >>> trig_finder.cd('170817529')
    >>> trig_finder
    <TriggerFtp: 170817529>
    >>> trig_finder.num_files
    128

Of course, you don't want to just list the files in a directory, you want to 
download them. Let's download all the catalog files for GRB 170817A:

    >>> trig_finder.get_cat_files('./', verbose=True)
    glg_bcat_all_bn170817529_v01.fit [==============================] 100.00%
    glg_scat_all_bn170817529_flnc_band_v00.fit [==============================] 100.00%
    glg_scat_all_bn170817529_flnc_comp_v00.fit [==============================] 100.00%
    glg_scat_all_bn170817529_flnc_plaw_v00.fit [==============================] 100.00%
    glg_scat_all_bn170817529_flnc_sbpl_v00.fit [==============================] 100.00%
    glg_scat_all_bn170817529_pflx_band_v00.fit [==============================] 100.00%
    glg_scat_all_bn170817529_pflx_comp_v00.fit [==============================] 100.00%
    glg_scat_all_bn170817529_pflx_plaw_v00.fit [==============================] 100.00%
    glg_scat_all_bn170817529_pflx_sbpl_v00.fit [==============================] 100.00%
    glg_tcat_all_bn170817529_v03.fit [==============================] 100.00%


Finding Continuous GBM Data
===========================
Now we want some continuous data. There aren't any trigger numbers for 
continuous data. Continuous CTIME and CSPEC are available in files that cover 
a whole day (in UTC) and TTE are offered in hourly files. To find the data you 
need, instead of a trigger number, you need to create a |ContinuousFtp|
object by specifying a time using Astropy Time:

    >>> from gdt.missions.fermi.gbm.finders import ContinuousFtp
    >>> from gdt.missions.fermi.time import Time
    >>> time = Time(587683338.0, format='fermi')
    >>> cont_finder = ContinuousFtp(time)
    >>> cont_finder
    <ContinuousFtp: 587683338.0>
    >>> cont_finder.num_files
    379

That's a whole lotta files in this directory. Most of them are TTE; remember 
that each hour has a TTE file (since the end of 2012) for each detector. Let's 
just list the CTIME that's available:

    >>> cont_finder.ls_ctime()
    ['glg_ctime_b0_190816_v00.pha',
     'glg_ctime_b1_190816_v00.pha',
     'glg_ctime_n0_190816_v00.pha',
     'glg_ctime_n1_190816_v00.pha',
     'glg_ctime_n2_190816_v00.pha',
     'glg_ctime_n3_190816_v00.pha',
     'glg_ctime_n4_190816_v00.pha',
     'glg_ctime_n5_190816_v00.pha',
     'glg_ctime_n6_190816_v00.pha',
     'glg_ctime_n7_190816_v00.pha',
     'glg_ctime_n8_190816_v00.pha',
     'glg_ctime_n9_190816_v00.pha',
     'glg_ctime_na_190816_v00.pha',
     'glg_ctime_nb_190816_v00.pha']

Now let's list the available TTE for this time. This will only list the TTE 
files in the directory that cover the relevant time:

    >>> cont_finder.ls_tte()
    ['glg_tte_b0_190816_21z_v00.fit.gz',
     'glg_tte_b1_190816_21z_v00.fit.gz',
     'glg_tte_n0_190816_21z_v00.fit.gz',
     'glg_tte_n1_190816_21z_v00.fit.gz',
     'glg_tte_n2_190816_21z_v00.fit.gz',
     'glg_tte_n3_190816_21z_v00.fit.gz',
     'glg_tte_n4_190816_21z_v00.fit.gz',
     'glg_tte_n5_190816_21z_v00.fit.gz',
     'glg_tte_n6_190816_21z_v00.fit.gz',
     'glg_tte_n7_190816_21z_v00.fit.gz',
     'glg_tte_n8_190816_21z_v00.fit.gz',
     'glg_tte_n9_190816_21z_v00.fit.gz',
     'glg_tte_na_190816_21z_v00.fit.gz',
     'glg_tte_nb_190816_21z_v00.fit.gz']

Similar to the trigger finder, you can use the same object to search at 
different times:

    >>> new_time = Time('2017-08-17T12:41:06.47', format='isot', scale='utc')

Now how about downloading the position history file for this time:

    >>> cont_finder.get_poshist('./', verbose=True)
    glg_poshist_all_170817_v01.fit [==============================] 100.00%


See :ref:`The FtpFinder Class` for more details on using data finders.

Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.finders
   :inherited-members:


