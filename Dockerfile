FROM regfe89/apache-ubuntu-18:8

COPY ./app /app

COPY modify.sh /app

RUN cd /app && ./modify.sh

RUN apt update && apt install make csh tcsh libx11-dev gcc-4.8 gfortran-4.8 gmt gmt-dcw gmt-gshhg gcc-multilib gfortran-multilib -y

RUN cd /usr/bin && ln -sf gcc-4.8 gcc && ln -sf gfortran-4.8 gfortran

RUN cd /app && ./install_software

RUN cd /app/tables && rm -f otl_FES2004.grid

COPY ./otl_FES2004.grid /app/tables/otl_FES2004.grid

ENV PATH="/app/com/:/app/gamit/bin:/app/kf/bin:${PATH}"
