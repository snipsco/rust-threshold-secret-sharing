export NDK_HOME=$HOME/Library/Android/sdk/ndk-bundle
export ANDROID_SDK_HOME=$HOME/Library/Android/sdk

cargo apk build

#rustc --target=arm-linux-androideabi \
#    -C linker=$ANDROID_SDK_HOME/ndk-bundle/toolchains/arm-linux-androideabi-4.9/prebuilt/darwin-x86_64/bin/arm-linux-androideabi-ld \
#    -C link-args="-pie -fPIE" src/lib.rs

