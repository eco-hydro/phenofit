Running in docker: 
```bash
docker run -it -p 8787:8787 `
    -v E:/ubuntu:/ubuntu `
    -v E:/ubuntu/site-library:/usr/local/lib/R/site-library `
    -v E:/ubuntu/etc:/usr/local/lib/R/etc `
    -v E:/github:/github `
    --name R-3.5.3 kongdd/phenofit bash 

# powershell docker rm @(docker ps -aq)
# powershell docker rmi @(docker images -f "dangling=true" -q)
```
