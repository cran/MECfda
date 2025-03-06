
# **Changelog**

## **Version 0.2.0**  

**Maintainer:** Heyang Ji  

**Previous Version:** 0.1.0  

### **New Features**  
- Added classes and functions required for basis expansion using **Functional Principal Component (FPC)** basis, including:  
  - **s4 classes**: `numeric_basis`, `numericBasis_series`.  
  - **s4 class methods**: `plot.numeric_series`.  
  - **Functions**: `FPC_basis_expansion`(), `numeric_basis_expansion()`, `numericBasisSeries2fun()`.  

### **Enhancements & Compatibility Improvements**  
- Enhanced existing functions to support the newly introduced FPC basis expansion:  
  - Added an option to use the FPC basis in **scalar-on-function regression functions**, `fcRegression` and `fcQR`.  
  - Updated **S3 methods** `predict.fcRegression` and `predict.fcQR` to support FPC basis.  
  - Extended the **generic function** `basis2fun` with new methods.  

### **New Utility Functions**  
- Introduced two functions for generating synthetic test datasets:  
  - `MECfda_simDataGen_fcReg`  
  - `MECfda_simDataGen_ME`  

### **Bug Fixes**  
- Fixed an issue in `bspline_basis_expansion` and `fourier_basis_expansion` where the numerical integration step did not properly handle **non-default integration intervals**.  

### **Other Improvements**  
- Added the function `MEM_X_hat()`, which returns $\hat X(t)$ used in `ME.fcRegression_MEM`.  
- Removed the package dependency on **stringr**.  
- Updated and refined the **help documentation** and **vignette** content.  
- Corrected **author roles** in the author list.  
