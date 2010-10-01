********       Include for common block trans  **********************

	Real P(Nlevel)
	Real T(Nlevel,Natm),Q(Nlevel,Natm)
        Real Pw(Nlevel,Natm),O3(Nlevel,Natm)
	Real DryTrans(0:Nlevel,Nchan,Nangle,Natm)
	Real WetTrans(0:Nlevel,Nchan,Nangle,Natm)
	Real  O3Trans(0:Nlevel,Nchan,Nangle,Natm)
	Real Tb(Nchan,Nangle,Natm)

	Common/trans/P,T,Q,Pw,O3,DryTrans,WetTrans,O3Trans,Tb

********       End Include for common block trans ********************
