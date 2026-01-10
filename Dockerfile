# Use Ubuntu 22.04 base
FROM ubuntu:22.04

# 2. Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
    PATH="/opt/conda/bin:$PATH"

# 3. System basics + libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    libnsl2 \
    git \
    libgl1-mesa-glx \
    libglu1-mesa \
    ca-certificates \
    build-essential \
    faketime \
    && rm -rf /var/lib/apt/lists/*

# 4. Install Miniforge (which includes Mamba)
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh && \
    bash Miniforge3-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniforge3-Linux-x86_64.sh

# 5. Copy environment.yml and build the environment
COPY environment.yml /project/environment.yml
WORKDIR /project
RUN mamba env create -f environment.yml && mamba clean -afy

# 6. Install FoldX
COPY bin/foldx /usr/local/bin/foldx_bin
RUN chmod +x /usr/local/bin/foldx_bin

# Create wrapper script that runs FoldX with faketime command
RUN echo '#!/bin/bash\nexec faketime "-1y" /usr/local/bin/foldx_bin "$@"' > /usr/local/bin/foldx && \
    chmod +x /usr/local/bin/foldx

ENV FAKETIME="-1y"

ENTRYPOINT ["/bin/bash", "-c"]