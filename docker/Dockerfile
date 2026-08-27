# pull official base image
FROM gcc:11-bullseye


# Create code dir
RUN mkdir /usr/src/app
WORKDIR /usr/src/app
RUN mkdir /usr/src/app/deps

# install dependencies
COPY --chmod=777 ./install/matio.sh ./deps/
RUN ./deps/matio.sh
COPY --chmod=777 ./install/hdf5.sh ./install/hdf5-1.14.2.tar.gz ./deps/
RUN ./deps/hdf5.sh
COPY --chmod=777 ./install/boost.sh ./deps/
RUN ./deps/boost.sh
# Bullseye EOL: point to archive, drop security repo (not mirrored there)
RUN sed -i 's/deb.debian.org/archive.debian.org/g' /etc/apt/sources.list \
    && sed -i '/bullseye-security/d' /etc/apt/sources.list \
    && echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99no-check-valid-until
RUN apt-get update && apt-get install libcgal-dev -y
COPY --chmod=777 ./install/cereal.sh ./deps/
RUN ./deps/cereal.sh
COPY --chmod=777 ./install/eigen.sh ./deps/
RUN ./deps/eigen.sh
# Fix for gfortran not found when installing openmpi
RUN ln -s /usr/bin/gfortran-10 /usr/bin/gfortran
RUN apt-get install -y libopenmpi-dev openmpi-bin

RUN apt install openssh-server -y
# Create a user “sshuser” and group “sshgroup”
ARG APP_UID=10295
RUN groupadd -g ${APP_UID} sshgroup && useradd -ms /bin/bash -g sshgroup -u ${APP_UID} sshuser
# Create sshuser directory in home
RUN mkdir -p /home/sshuser/.ssh
# Copy the ssh public key in the authorized_keys file. The idkey.pub below is a public key file you get from ssh-keygen. They are under ~/.ssh directory by default.
COPY floe_rsa.pub /home/sshuser/.ssh/authorized_keys
# change ownership of the key file. 
RUN chown sshuser:sshgroup /home/sshuser/.ssh/authorized_keys && chmod 600 /home/sshuser/.ssh/authorized_keys
# Start SSH service
RUN service ssh start

# Install visu dependencies
RUN apt-get install python3-pip -y
COPY requirements.txt .
RUN pip install -r requirements.txt
RUN apt-get update \
 && apt-get install -y --no-install-recommends ffmpeg \
 && rm -rf /var/lib/apt/lists/*

# entrypoint
COPY --chmod=777 ./entrypoint.sh /usr/src/app/entrypoint.sh
WORKDIR /usr/src/app/floedyn
ENTRYPOINT ["/usr/src/app/entrypoint.sh"]
