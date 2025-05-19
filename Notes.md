# Data Checks

- Abramovich data (A0244002) are cut in AZURE2 files: the real ones start at 0.05 MeV, the AZURE2 at 0.12 MeV. The reason is the screening effect.

- Barnard data (A1269002) are not usedd because the newer Paneru data are available.

- Elwyn differential (F0012002) had lower uncertainties w.r.t. EXFOR.

- Elwyn differential (F0012003) had much lower uncertainties w.r.t. EXFOR.

- Elwyn integrated data (F0012005) has a small shift down of 5 keV with respect to EXFOR and higher uncertainties for some reason. We fixed it and are using the new EXFOR data since the difference is minimal and the correction might have been done for unknow reason by Ian.

- Fasoli data (D0135002) had smaller uncertainties w.r.t. EXFOR and was shifted 30 keV down for some reason. We use the EXFOR now.

- Harrison data (C1003002) is the same as the EXFOR one.

- Ivanovic data (A1014010) used 3He beam instead of 4He one so the energies were converted accordingly.

- McCray data (A1410002, A1410003 and A1410004) are problematic since Ian seems to have many more data with respect to these available on EXFOR. However, EXFOR data points are more dense, the Ian's one instead seems to have more angles instead. The A1410004 seems problematic to fit: James provided me his digitized data and they seem to be much better.

- Schenk data (F0049002) are okay.

- Spiger data (A1094004) was ok.

- Spiger data (A1094007) are in lab angles on EXFOR respect to CM in Ian's data. Additionally, the A1094006 dataset is available but some of the data at low energies are too low with respect to the AZURE2 prediction thus they are ignored.

- Tumino data (O1221002) is the same as the EXFOR, but they had to be converted to lab energies for AZURE2.

- Sometimes the errors in EXFOR are null, thus we assume a conservative 10% for these.

