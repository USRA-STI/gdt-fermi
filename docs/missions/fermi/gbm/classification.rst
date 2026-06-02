.. _gbm-classification:
.. |Trigdat| replace:: :class:`~gdt.missions.fermi.gbm.trigdat.Trigdat`
.. |Classification| replace:: :class:`~gdt.missions.fermi.gbm.classification.Classification`
.. |classify()| replace:: :meth:`~gdt.missions.fermi.gbm.classification.Classification.classify`
.. |classify_trigdat()| replace:: :meth:`~gdt.missions.fermi.gbm.classification.Classification.classify_trigdat`

*********************************************************************************
Fermi GBM Classification Algorithm (:mod:`gdt.missions.fermi.gbm.classification`)
*********************************************************************************

The Fermi GBM has an on-board classifier for high-energy transients detected by
the flight software.  This classifier is operated in real-time on the spacecraft
and it uses information such as the on-board localization and hardness ratio to
determine if the trigger is astrophysical, and if so, what is the likely source
of the trigger.  

.. Note::
  This GDT implementation of the GBM classification algorithm is
  a close approximation of the algorithm operating on the spacecraft; however 
  there may be small differences compared to the exact flight software algorithm.


While there are a couple of ways to employ this algorithm, we will use the 
method that reads the required information from a trigdat (trigger data) file
(see :ref:`this documentation <gbm-trigdat>` for more info on working
with trigdat). The alternative method requires cobbling together the requisite
inputs, which an advanced user may wish to do when finding a potential signal
in continuous GBM that did not produce an on-board trigger (see the |classify()|
method for more info about the inputs).


To start, we will create and initialize our |Classification| object:

    >>> from gdt.missions.fermi.gbm.classification import Classification
    >>> classifier = Classification()


Next, we open a trigdat file:

    >>> from gdt.core import data_path
    >>> from gdt.missions.fermi.gbm.trigdat import Trigdat
    >>> filepath = data_path.joinpath('fermi-gbm/glg_trigdat_all_bn170101116_v01.fit')
    >>> trigdat = Trigdat.open(filepath)
    
Now we can pass the |Trigdat| object to the |Classification| object, which 
performs an  approximation of the source classification routine used on-board 
GBM.  During this step we apply the following settings:

* ``use_loc=True`` performs the classification using the detectors that are 
  closest to the source location. When ``False`` the routine will instead sort 
  the detectors by the largest count rates in the 50-300 keV energy band.

* ``verbose=True`` prints the trigdat information used as part of the 
  classification algorithm. Setting this to ``False``` or ignoring this argument 
  will silence this output.

    >>> det_nums = classifier.classify_trigdat(trigdat, use_loc=True, verbose=True)
    Trigger Algorithm = 6
    Location = (RA 59.167, Dec 9.333) deg
    Location Err = 5.317 deg
    Rates
     [[ 64. 279. 223. 246. 269.  46.  38.  72.]
     [ 54. 283. 219. 228. 235.  40.  43.  57.]
     [ 42. 202. 197. 194. 202.  40.  56.  37.]
     [ 67. 290. 218. 191. 184.  36.  36.  74.]
     [ 43. 216. 202. 204. 176.  39.  67.  47.]
     [ 35. 154. 183. 198. 205.  45.  90.  29.]
     [ 82. 364. 269. 243. 265.  48.  65.  20.]
     [ 99. 410. 291. 245. 253.  45.  44.  46.]
     [ 70. 392. 236. 211. 195.  39.  62.  36.]
     [ 82. 518. 379. 375. 342.  51.  61.  16.]
     [ 91. 483. 396. 371. 348.  46.  68.  22.]
     [ 74. 438. 303. 307. 279.  49.  32.  47.]
     [623. 399. 525. 224.  36.  30.  21.  95.]
     [702. 423. 423. 172.  28.  23.  21.  68.]]
    Background
     [[ 63. 276. 195. 145. 153.  39.  33.  65.]
     [ 56. 277. 205. 148. 147.  40.  43.  54.]
     [ 45. 206. 177. 145. 156.  38.  64.  41.]
     [ 71. 303. 210. 160. 159.  42.  40.  82.]
     [ 45. 207. 181. 149. 160.  40.  76.  38.]
     [ 38. 145. 155. 138. 172.  47.  98.  27.]
     [ 72. 325. 201. 157. 141.  34.  69.  21.]
     [ 88. 387. 240. 161. 140.  32.  45.  45.]
     [ 68. 372. 232. 163. 141.  34.  61.  33.]
     [ 70. 403. 253. 167. 137.  31.  61.  16.]
     [ 86. 382. 257. 169. 136.  31.  59.  19.]
     [ 63. 349. 227. 157. 143.  35.  31.  52.]
     [550. 402. 488. 240.  40.  29.  22.  91.]
     [412. 347. 387. 168.  31.  24.  18.  75.]]


The |classify_trigdat()| method returns a list with the detector numbers for 
the top 2 detectors closest to the source location. Let's check which ones 
they are:

    >>> from gdt.missions.fermi.gbm.detectors import GbmDetectors
    >>> print("\nClosest detectors:")
    >>> for num in det_nums:
    >>>     print("  {}".format(GbmDetectors.from_num(num).name))    
    Closest detectors:
      n9
      na

We can also inspect the top two source classifications, ranked according to 
probability:

    >>> print("\nTop 2 classifications (this calculation):")
    >>> rank = classifier.rankings()[:2]
    >>> for r in classifier.rankings()[:2]:
    >>>     print("  prob = {:.2f} {}".format(r[1], r[0].name))   
    Top 2 classifications (this calculation):
      prob = 0.86 GRB
      prob = 0.10 GROJ_0422_32


For comparison, let's also show the classifications determined by the flight 
software on-board the spacecraft:

    >>> fsw_loc = trigdat.fsw_locations[-1]
    >>> print("\nTop 2 classifications (flight):")
    >>> for c, r in [fsw_loc.top_classification, fsw_loc.next_classification]:
    >>>     print("  prob = {:.2f} {}".format(r, c))    
    Top 2 classifications (flight):
      prob = 0.86 GRB
      prob = 0.09 GROJ422

And since we are curious, we can see the full breakdown of classification 
probabilities:

    >>> print(classifier)
    Classification(
                   ERROR               0.0000
                   UNRELIABLE_LOCATION 0.0000
                   LOCAL_PARTICLES     0.0000
                   BELOW_HORIZON       0.0000
                   GRB                 0.8633
                   GENERIC_SGR         0.0044
                   GENERIC_TRANSIENT   0.0359
                   DISTANT_PARTICLES   0.0000
                   SOLAR_FLARE         0.0000
                   CYG_X1              0.0000
                   SGR_1806_20         0.0000
                   GROJ_0422_32        0.0964
                   RESERVED_1          0.0000
                   RESERVED_2          0.0000
                   RESERVED_3          0.0000
                   RESERVED_4          0.0000
                   RESERVED_5          0.0000
                   RESERVED_6          0.0000
                   RESERVED_7          0.0000
                   TGF                 0.0000)

    
Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.classification

   :inherited-members:



