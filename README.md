# The ```rcdkTools``` package

This package provides some functions to simplify the fingerprint calculation for
larger sets of molecules. It is mainly indeted to be a wrapper around the
[```rcdk-package```](https://github.com/rajarshi/cdkr).

## TODO:
- Add convenient functions to handle counting fps
  - store counting fps matrices
  - store counting fps hashed lists
- Add tools to handle molecular descriptors
- Add parallel calculation of fingerprints and descriptors

## My modified ```rcdk-package```

Currently not support by the original package:
- **substructure** fingerprint provided by [CDK](https://github.com/cdk/cdk)
  - In this way the MACCS counting fingerprints can be used. 
  - For that the corresponding SMARTS pattern must be [downloaded](https://github.com/cdk/cdk/blob/4004eb64fd7e94a0da674ae2c0eedba79fda425f/descriptor/fingerprint/src/main/resources/org/openscience/cdk/fingerprint/data/SMARTS_countable_MACCS_keys.txt).
- additional **circular** fingerprint type as supported by CDK
  - The original ```rcdk``` only supports *ECFP6*
  - My extension provides access to the [all types provided by CDK](https://github.com/cdk/cdk/blob/4004eb64fd7e94a0da674ae2c0eedba79fda425f/descriptor/fingerprint/src/main/java/org/openscience/cdk/fingerprint/CircularFingerprinter.java#L99)

### Installation

Open a R terminal:

```R
library(devtools)
install_github("bachi55/cdkr", subdir="rcdklibs", ref="support_cdk_substructure_fps")
install_github("bachi55/cdkr", subdir="rcdk", ref="support_cdk_substructure_fps")
```

Please note: This overwrites the possible already existing installations.
