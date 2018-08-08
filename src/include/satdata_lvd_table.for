      character*6 satellite_ids(maxsat)
      character*3 satellite_types(maxtype,maxsat)
      character*3 satellite_channels(maxtype,maxchannel,maxsat)

      data satellite_ids/
     1  'goes08',
     2  'meteos',
     3  'goes10',
     4  'gmssat',
     5  'goes12',
     6  'goes09',
     7  'goes11',
     8  'noaapo',
     9  'mtsat',
     .  'fy',
     1  'goeswe',
     2  'goesea',
     3  'goesne',
     4  'goesnw',
     5  'goesse',
     6  'goessw',
     7  'goes16',
     8  'goes17',
     9  'coms',
     .  'him8',
     1  'him9'/

      data satellite_types/
     1              'gvr','wfo','cdf','rll',
     2              '   ','   ','   ','rll',
     3              'gvr','wfo','   ','   ',
     4              '   ','hko','twn','rll',
     5              'gvr','wfo','cdf','rll',
     6              'gvr','   ','cdf','   ',
     7              'gvr','wfo','cdf','rll',
     8              'ncp','   ','   ','rll',
     9              '   ','   ','   ','rll',
     .              '   ','   ','   ','rll',
     1              'gvr','wfo','cdf','rll',
     2              'gvr','wfo','cdf','rll',
     3              'gvr','wfo','cdf','rll',
     4              'gvr','wfo','cdf','rll',
     5              'gvr','wfo','cdf','rll',
     6              'gvr','wfo','cdf','rll',
     7              'gnp','   ','   ','   ',
     8              'gnp','   ','   ','   ',
     9              'cms','   ','   ','   ',
     .              'rll','nll','jma','   ',
     1              'rll','nll','jma','   '/

      data satellite_channels/
     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              'vis','4u ','wv ','11u','12u','   ',  !end for goes08

     1              '   ','   ','   ','   ','   ','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              'vis','4u ','wv ','10p','12u','   ',  !end for meteos

     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goes10

     1              '   ','   ','   ','   ','   ','   ',
     2              'vis','   ','wv ','11u','12u','   ',
     3              'vis','   ','wv ','11u','12u','   ',
     4              'vis','   ','wvp','11u','12u','   ',  !end for gmssat
 
     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goes12

     1              'vis','4u ','wv ','11u','12u','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goes09

     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goes11

     1              'vis','i39','   ','i11','   ','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              'vis','4u ','   ','11u','   ','   ',  !end for noaapo

     1              '   ','   ','   ','   ','   ','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              '   ','   ','   ','10p','   ','   ',  !end for mtsat

     1              '   ','   ','   ','   ','   ','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              '   ','   ','   ','10p','   ','   ',  !end for fy

     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goeswe

     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goesea

     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goesne

     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goesnw

     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goesse

     1              'vis','4u ','wv ','11u','12u','   ',
     2              'vis','i39','iwv','i11','i12','   ',
     3              'vis','4u ','wv ','11u','12u','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goessw

     1              'vis','4u ','wv ','11u','12u','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goes16

     1              'vis','4u ','wv ','11u','12u','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for goes17

     1              'vis','4u ','wv ','11u','12u','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for coms

     1              'vis','4u ','wv ','10p','12u','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              '   ','   ','   ','   ','   ','   ',  !end for him8

     1              'vis','4u ','wv ','10p','12u','   ',
     2              '   ','   ','   ','   ','   ','   ',
     3              '   ','   ','   ','   ','   ','   ',
     4              '   ','   ','   ','   ','   ','   '/  !end for him9
