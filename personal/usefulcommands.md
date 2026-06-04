# remove ^M in files for nvim
```bash
:%s/
//g
```

# pipe to clipboard
## just clipboard
```bash
cat hello | clip.exe
```
## output to terminal and clipboard
```bash
cat hello | tee clip.exe
```