class ProfilesController < ApplicationController

  def new
    @profile = Profile.new
    @profiles = Profile.all
  end

  def create

    # save image file
    uploaded_io = params["profile"]["profile_image"]
    File.open(Rails.root.join('public','profiles',
      uploaded_io.original_filename), 'wb') do |file|
        file.write(uploaded_io.read)
      end

    # grab params for profile
    profile_params = params["profile"]
    name = profile_params["name"]
    email = profile_params["email"]
    profile_image = uploaded_io.original_filename
    description = profile_params["description"]
    @profile = Profile.new(name:name,
                          email:email,
                          description: description,
                          img_url:profile_image)
    @profile.save
    # binding.pry
    redirect_to maintenance_index_path
  end

  def edit

  end

  def updated

  end

  def destroy
    profile = Profile.find(params[:id])
    profile.destroy
    redirect_to new_profile_path
  end

  end
