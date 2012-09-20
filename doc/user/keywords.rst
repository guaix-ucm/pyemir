
FITS Keywords
=============
The description of the keywords follows a convention found in other FITS 
keyword dictionaries, for example the list in http://fits.gsfc.nasa.gov/fits_dictionary.html. 
The keyword name is expressed, with the reference to the paper where it is 
included. Following the type of HDU where the keyword can appear. The value 
shows the kind of variable represented by the keyword. The comment is a 
example of the comment associated with the keyword and definition is a 
explanation in detail of the usage of the keyword.

Primary header
--------------
==== ========  ========================  ============================================
Type Keyword   Example                    Explanation
==== ========  ========================  ============================================
 L   SIMPLE     T                        Standard FITS format
 I   BITPIX     16                       One of -64,-32,8,16,32
 I   NAXIS      2                        # of axes in frame
 I   NAXIS1    2048                      # of pixels per row
 I   NAXIS2    2048                      # of rows
 S   ORIGIN    'GTC'                     FITS originator
 S   OBSERVAT  'ORM'                     Observatory
 S   TELESCOP  'GTC'                     The telescope
 S   INSTRUME  'EMIR'                    The instrument
 S   OBJECT    'NGC 4594'                Target designation
 S   OBSERVER  'OBSERVER'                Who adquired the data
 S   DATE-OBS  '2012-09-20T12:00:11.50'  Date of the start of the observation
 S   DATE      '2012-09-20T12:14:12.78'  Date the file was written
==== ========  ========================  ============================================

Required by the pipeline
------------------------

==== ========  ========================  ============================================
Type Keyword   Example                    Explanation
==== ========  ========================  ============================================
 R   AIRMASS   1.1908                    Mean airmass of the observation
 R   MJD-OBS   72343.34324               Modified JD of the start of the observation
 S   IMAGETYP  'FLAT'                    Type of the image
 S   OBSTYPE   'FLATON'                  Type of observation
 L   READPROC   T                        The frame has been preprocessed after readout
 S   READMODE  'RAMP'                    The readmode used to adquire the image
 I   READSAMP   12                       Number of samples taken (equal to NAXIS1 if NAXIS is 3)
 R   EXPOSED                             Photometric time?
 R   DARKTIME
 R   EXPTIME
 R   ELAPSED                             Time between resets?
 I   OBSID      567                      Identifier of the observing block
 S   OBSMODE   'DITHER_IMAGES'           Identifier of the observing mode
 I   OBSEXPN    30                       # of exposures during the block
 I   OBSSEQN    12                       # of the image in the sequence
==== ========  ========================  ============================================

Coordinate system
-----------------

Checksum convention
-------------------
The CHECKSUM and DATASUM keywords that are embedded in the FITS header are used to verify the integrity of the HDU.

See http://fits.gsfc.nasa.gov/registry/checksum.html

==== ========  ========================  ============================================
Type Keyword   Example                    Explanation
==== ========  ========================  ============================================
 S   CHECKSUM  'ADFASASDLIEXV'           HDU checksum 
 S   DATASUM   '1929302939392            data unit checksum
==== ========  ========================  ============================================
