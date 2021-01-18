from frospy.spectrum.app import spectrum

data = '/quanta1/home/simons/splitting/modes/00s13/self_coupling/R-Comp/12_plus/syn/010196C.ahx.syn'
syn = '/quanta1/home/simons/splitting/modes/00t14/self_coupling/R-Comp/12_plus/syn/010196C.ahx.syn'

spectrum(data=data, syn=syn, fw=[1,2], tw=[3,60], minispec=True
