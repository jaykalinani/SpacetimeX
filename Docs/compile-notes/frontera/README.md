# Install CarpetX with Spack (Frontera)

* Use interactive session

    - Compile CPU version: `idev -m 120`

    - Compile GPU version: `idev -p rtx-dev -m 120`

* Download spack

    - `git clone -c feature.manyFiles=true https://github.com/spack/spack.git`
    
    - `git checkout relesases/v0.20`

* Download CarpetX and SpacetimeX

    ```bash
    curl -kLO https://raw.githubusercontent.com/gridaphobe/CRL/master/GetComponents
    chmod a+x GetComponents
    ./GetComponents --root Cactus --parallel --no-shallow https://raw.githubusercontent.com/lwJi/SpacetimeX/main/Docs/thornlist/spacetimex.th
    ```

## Intel-Oneapi version (`oneapi@2023.1.0`)

* Load `intel/23.1.0`

    - `module load intel/23.1.0`

* Setup spack

    - `. share/spack/setup-env.sh`
    
    - `spack compiler find`

* Create a dir where you want put `view` in (say `/work2/.../username/frontera/SpackView/oneapi`)

* Replace the last line of `oneapi-23.1.0/spack.yaml` with your own dir (say `/work2/.../username/frontera/SpackView/oneapi/view`)

* Replace the dir `/work2/08708/liwei/frontera/SpackView/oneapi/view` in `configs/config_frontera_oneapi.cfg` (with say `/work2/.../username/frontera/SpackView/oneapi/view`)

* Install other required packages

    ```bash
    env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 compiler find
    env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 concretize --force
    env TMPDIR=$WORK/tmp spack --env-dir ./oneapi-23.1.0 install --fail-fast
    ```

* Install SpacetimeX

    ```bash
    cd Cactus
    gmake SpacetimeX-oneapi options=repos/SpacetimeX/Docs/compile-notes/frontera/configs/config_frontera_oneapi.cfg
    cp repos/SpacetimeX/Docs/thornlist/spactimex.th configs/SpacetimeX-oneapi/ThornList
    gmake -j24 SpacetimeX-oneapi
    ```


## CUDA version (`cuda@11.8.0`)

* Make sure rerun `spack install gcc@11.2.0 %gcc@4.8.5` again on `rtx-dev` or `rtx`

* Create a dir where you want put `view` in (say `/work2/.../username/frontera/SpackView/cuda`)

* Replace the last line of `cuda-11.8.0/spack_yaml` with your dir (say `/work2/.../username/frontera/SpackView/cuda/view`)

* Replace the dir `/work2/08708/liwei/frontera/SpackView/cuda/view` in `configs/config_frontera_cuda.cfg` (with say `/work2/.../username/frontera/SpackView/cuda/view`)

* Install other required packages

    ```bash
    env TMPDIR=$WORK/tmp spack --env-dir ./cuda-11.8.0 compiler find
    env TMPDIR=$WORK/tmp spack --env-dir ./cuda-11.8.0 concretize --force
    env TMPDIR=$WORK/tmp spack --env-dir ./cuda-11.8.0 install --fail-fast
    ```

* Install SpacetimeX

    ```bash
    spack load gcc@11.2.0
    spack load cuda@11.8.0
    cd Cactus
    gmake SpacetimeX-cuda options=repos/SpacetimeX/Docs/compile-notes/frontera/configs/config_frontera_cuda.cfg
    cp repos/SpacetimeX/Docs/thornlist/spactimex.th configs/SpacetimeX-cuda/ThornList
    gmake -j16 SpacetimeX-cuda
    ```
