# MODELLER 10.8 Installation Notes

Installation notes for MODELLER 10.8 on Ubuntu (Linux x86_64) using Miniconda.

- **Platform:** Ubuntu, 64-bit (amd64)
- **Python:** 3.13 (conda base environment)
- **Method used:** conda (recommended)

---

## Method 1: conda (recommended)

### 1. Add the Salilab channel

```bash
conda config --add channels salilab
```

### 2. Install MODELLER with the license key

```bash
KEY_MODELLER=MODELIRANJE conda install modeller
```

Setting `KEY_MODELLER` before `conda install` writes the license key into MODELLER's config automatically. Answer `y` when prompted to proceed.

> Replace `MODELIRANJE` with your own key if you have a different one. Academic keys are free; request one from the [MODELLER registration page](https://salilab.org/modeller/registration.html).

### 3. Verify the import

```bash
python -c "import modeller; print(modeller.__file__)"
```

Expected output (path will vary with your username and conda location):

```
/home/<user>/miniconda3/lib/python3.13/site-packages/modeller/__init__.py
```

### 4. Verify the license

```bash
python -c "from modeller import *; env = Environ(); print('Modeller OK')"
```

If it reports a license error, edit `.../site-packages/modeller/config.py` and set:

```python
license = 'MODELIRANJE'
```

---

## Running MODELLER

MODELLER is a Python library, so there is **no `modeller` command**. Write a script and run it inside the conda environment:

```bash
python foo.py
```

For a first script, see the [MODELLER tutorial](https://salilab.org/modeller/tutorial/). The basic comparative modelling workflow is `align2d` followed by `model-single`.

---

## Method 2: `.deb` package (alternative)

```bash
cd ~/Downloads
wget https://salilab.org/modeller/10.8/modeller_10.8-1_amd64.deb
ls -l modeller_10.8-1_amd64.deb      # should be tens of MB, not 0 bytes
sudo env KEY_MODELLER=MODELIRANJE dpkg -i modeller_10.8-1_amd64.deb
```

- Installs to `/usr/lib/modeller10.8/`.
- Run scripts with the system interpreter (`python3 foo.py`) or the bundled launcher (`mod10.8 foo.py`).
- It is **not** importable from conda Python, because it uses the system `python3`.
- Optional: verify the download with the GPG signature (`modeller_10.8-1_amd64.deb.sign`), after importing the MODELLER signing key:

  ```bash
  gpg --verify modeller_10.8-1_amd64.deb.sign modeller_10.8-1_amd64.deb
  ```

---

## Troubleshooting

| Problem | Cause | Fix |
|---|---|---|
| `dpkg-deb: error: unexpected end of file in archive magic version number` | The `.deb` download was empty or truncated | Check with `ls -l` and `file`. Delete the file and re-download, or use the conda method |
| Downloaded file is only a few hundred bytes | The `.deb.sign` link is a GPG signature, not the package | Download the URL without `.sign` |
| `modeller: command not found` | MODELLER is a Python library, not a standalone command | Use `python foo.py` |
| `import modeller` works but license error appears | Key not written to config | Set `license = '...'` in `site-packages/modeller/config.py` |
| Anaconda notice about Python 3.10 builds ending | General channel warning | Not related to this install; safe to ignore on Python 3.11+ |

---

## References

- [MODELLER 10.8 release notes](https://salilab.org/modeller/release.html)
- [MODELLER tutorial](https://salilab.org/modeller/tutorial/)
- [MODELLER registration (license key)](https://salilab.org/modeller/registration.html)
