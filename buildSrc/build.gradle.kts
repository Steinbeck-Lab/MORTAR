plugins {
    // This version of Gradle expects version "4.2.1" of the `kotlin-dsl`
    id("org.gradle.kotlin.kotlin-dsl") version "4.2.1"
}

repositories {
    // Needed to resolve the Kotlin-DSL plugin
    gradlePluginPortal()
}
