c  requires lapsparms.for
        integer igrid,jgrid,kgrid,kdim,ksurf,imaps_nest,jmaps_nest
        integer cloud_levels,layers,datainit
        real mdf

c reference lapsparms.for parameters for grid dimensions to consolodate
c                                                               changes
        parameter (igrid = NX_L)
        parameter (jgrid = NY_L)
        parameter (kgrid = NZ_L)
        parameter (kdim = NZ_L)

c                               !used for sh levels in analyses
        parameter (ksurf = 1)
c               ! used for glst for surface data (one vert dim)
c        parameter (imaps_nest = 11)
c        parameter (jmaps_nest = 11)
        parameter (imaps_nest = igrid)
        parameter (jmaps_nest = jgrid)
        parameter (cloud_levels = 42    )
        parameter (layers=3 )
        parameter (datainit = igrid*jgrid*kdim)
        parameter (mdf = r_missing_data) ! missing data flag from lapsparms.for

        character*31 cloud_ext
        data cloud_ext /'lc3'/
