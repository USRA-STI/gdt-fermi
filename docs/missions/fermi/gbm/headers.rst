.. _gbm-headers:

**************************************************************
Fermi GBM FITS Headers (:mod:`gdt.missions.fermi.gbm.headers`)
**************************************************************
This module defines all of the FITS headers for the public data files. While
these classes are not usually directly called by the user, we may load one up
and see the contents and default values.  For example, here is the set of 
header definitions for triggered PHAII (CTIME and CSPEC) files:

    >>> from gdt.missions.fermi.gbm.headers import PhaiiTriggerHeaders
    >>> hdrs = PhaiiTriggerHeaders()
    >>> hdrs
    <PhaiiTriggerHeaders: 4 headers>
    
Here is the ``PRIMARY`` header and default values (retrieved by index):

    >>> hdrs[0]
    CREATOR = 'Gamma-ray Data Tools 2.0.0' / Software and version creating file     
    FILETYPE= '' / Name for this type of FITS file                                  
    FILE-VER= '1.0.0   '           / Version of the format for this filetype        
    TELESCOP= 'GLAST   '           / Name of mission/satellite                      
    INSTRUME= 'GBM     '           / Specific instrument used for observation       
    DETNAM  = '' / Individual detector name                                         
    OBSERVER= 'Meegan  '           / GLAST Burst Monitor P.I.                       
    ORIGIN  = 'GIOC    '           / Name of organization making file               
    DATE    = '2023-02-12T22:18:54.906' / file creation date (YYYY-MM-DDThh:mm:ss UT
    DATE-OBS= '' / Date of start of observation                                     
    DATE-END= '' / Date of end of observation                                       
    TIMESYS = 'TT      '           / Time system used in time keywords              
    TIMEUNIT= 's       '           / Time since MJDREF, used in TSTART and TSTOP    
    MJDREFI =                51910 / MJD of GLAST reference epoch, integer part     
    MJDREFF = '7.428703703703703e-4' / MJD of GLAST reference epoch, fractional part
    TSTART  =                  0.0 / [GLAST MET] Observation start time             
    TSTOP   =                  0.0 / [GLAST MET] Observation stop time              
    FILENAME= '' / Name of this file                                                
    DATATYPE= '' / GBM datatype used for this file                                  
    TRIGTIME=                  0.0 / Trigger time relative to MJDREF, double precisi
    OBJECT  = '' / Burst name in standard format, yymmddfff                         
    RADECSYS= 'FK5     '           / Stellar reference frame                        
    EQUINOX =               2000.0 / Equinox for RA and Dec                         
    RA_OBJ  =                  0.0 / Calculated RA of burst                         
    DEC_OBJ =                  0.0 / Calculated Dec of burst                        
    ERR_RAD =                  0.0 / Calculated Location Error Radius               


And here is the ``SPECTRUM`` header and default values:

    >>> hdrs['SPECTRUM']
    EXTNAME = 'SPECTRUM'           / name of this binary table extension            
    TELESCOP= 'GLAST   '           / Name of mission/satellite                      
    INSTRUME= 'GBM     '           / Specific instrument used for observation       
    DETNAM  = '' / Individual detector name                                         
    OBSERVER= 'Meegan  '           / GLAST Burst Monitor P.I.                       
    ORIGIN  = 'GIOC    '           / Name of organization making file               
    DATE    = '2023-02-12T22:18:54.906' / file creation date (YYYY-MM-DDThh:mm:ss UT
    DATE-OBS= '' / Date of start of observation                                     
    DATE-END= '' / Date of end of observation                                       
    TIMESYS = 'TT      '           / Time system used in time keywords              
    TIMEUNIT= 's       '           / Time since MJDREF, used in TSTART and TSTOP    
    MJDREFI =                51910 / MJD of GLAST reference epoch, integer part     
    MJDREFF = '7.428703703703703e-4' / MJD of GLAST reference epoch, fractional part
    TSTART  =                  0.0 / [GLAST MET] Observation start time             
    TSTOP   =                  0.0 / [GLAST MET] Observation stop time              
    TRIGTIME=                  0.0 / Trigger time relative to MJDREF, double precisi
    OBJECT  = '' / Burst name in standard format, yymmddfff                         
    RADECSYS= 'FK5     '           / Stellar reference frame                        
    EQUINOX =               2000.0 / Equinox for RA and Dec                         
    RA_OBJ  =                  0.0 / Calculated RA of burst                         
    DEC_OBJ =                  0.0 / Calculated Dec of burst                        
    ERR_RAD =                  0.0 / Calculated Location Error Radius               
    FILTER  = 'none    '           / The instrument filter in use (if any)          
    AREASCAL=                  1.0 / No special scaling of effective area by channel
    BACKFILE= 'none    '           / Name of corresponding background file (if any) 
    BACKSCAL=                  1.0 / No scaling of background                       
    CORRFILE= 'none    '           / Name of corresponding correction file (if any) 
    CORRSCAL=                  1.0 / Correction scaling file                        
    RESPFILE= 'none    '           / Name of corresponding RMF file (if any)        
    ANCRFILE= 'none    '           / Name of corresponding ARF file (if any)        
    SYS_ERR =                  0.0 / No systematic errors                           
    POISSERR=                    T / Assume Poisson Errors                          
    GROUPING=                    0 / No special grouping has been applied           
    HDUCLASS= 'OGIP    '           / Conforms to OGIP standard indicated in HDUCLAS1
    HDUCLAS1= 'SPECTRUM'           / PHA dataset (OGIP memo OGIP-92-007)            
    HDUCLAS2= 'TOTAL   '           / Indicates gross data (source + background)     
    HDUCLAS3= 'COUNT   '           / Indicates data stored as counts                
    HDUCLAS4= 'TYPEII  '           / Indicates PHA Type II file format              
    HDUVERS = '1.2.1   '           / Version of HDUCLAS1 format in use              
    CHANTYPE= 'PHA     '           / No corrections have been applied               
    DETCHANS=                    0 / Total number of channels in each rate          
    EXTVER  =                    1 / Version of this extension format               

See :external:ref:`Data File Headers<core-headers>` for more information about 
creating and using FITS headers.
    
Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.headers
   :inherited-members:


