FROM python:3.11-bookworm

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    gfortran \
    libopenblas-dev \
    libffi-dev \
    libssl-dev \
    libcairo2 \
    libcairo2-dev \
    libpango-1.0-0 \
    libpango1.0-dev \
    libgdk-pixbuf-2.0-0 \
    libgdk-pixbuf-2.0-dev \
    shared-mime-info \
    curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .

RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Issue #66: bake the building commit into the image so generated reports can
# record which build produced their numbers. Passed by the GitHub Actions
# workflow; defaults to "unknown" for local builds.
ARG GIT_SHA=unknown
ENV MOLAOP_IMAGE_SHA=$GIT_SHA

RUN mkdir -p /app/uploads /app/temp && \
    useradd --create-home --shell /bin/bash appuser && \
    chown -R appuser:appuser /app/uploads /app/temp

USER appuser

EXPOSE 5000

CMD ["python", "app.py"]
