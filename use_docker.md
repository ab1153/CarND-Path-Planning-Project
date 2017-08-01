
* **run the following command in powershell**
```
> Set-NetConnectionProfile -interfacealias "vEthernet (DockerNAT)" -NetworkCategory Private
```
* **make the drive sharable**


* **run the docker**
```
> docker run -it -p 4567:4567 -v d:\side\CarND-Path-Planning-Project:/work udacity/controls_kit /bin/bash
```