Pickled files created using cPickle.

Example - inside `GW_link_Integrated_build.py`:

```python
res = 5000
hgu, hu_raster_files_reproj = \
        campaspe_mesh.build_mesh_and_set_properties(tr_model,
                                                    hu_raster_path,
                                                    hgu_props,
                                                    resolution=int(res)
                                                    )

# for tests - safe to remove once done
print(type(hgu), type(hu_raster_files_reproj))
import cPickle as pickle
with open('hgu_dump.pkl', 'wb') as hgu_dump:
    pickle.dump(hgu, hgu_dump)

with open('hgu_raster_files_reproj_dump.pkl', 'wb') as rpj:
    pickle.dump(hu_raster_files_reproj, rpj)

raise RuntimeError("refactoring stop")

## End for tests
```
