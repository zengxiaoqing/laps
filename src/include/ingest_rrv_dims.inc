c These values specify the rrv file size.
c     integer   nxv01,nyv01,nzv01
c     parameter (nxv01=121,nyv01=121,nzv01=21)

c Presently there are really only 3 radars for rrv ingest.
      integer   max_radars_rrv
      parameter (max_radars_rrv = 4)
c Allow processing of up to 3 files per ingest_rrv run
      integer   mxf
      parameter (mxf=3)
