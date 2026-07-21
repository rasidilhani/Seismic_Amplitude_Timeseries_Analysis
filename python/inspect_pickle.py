from pathlib import Path
import pickle

p = Path(r'E:\data\seismic_amplitude_timeseries_out\2019\2019.001\WIZ.NZ\2019.001.WIZ.RSAM.p')
print('exists', p.exists())
print('size', p.stat().st_size if p.exists() else None)
if p.exists():
    with open(p, 'rb') as f:
        obj = pickle.load(f)
    print(type(obj))
    if hasattr(obj, 'shape'):
        print('shape', obj.shape)
    if hasattr(obj, '__len__'):
        print('len', len(obj))
    print('repr', repr(obj)[:1000])
